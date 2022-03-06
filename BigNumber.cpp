#include "BigNumber.h"
#include <cmath>
#include <climits>
#include <ranges>
#include <algorithm>
#include <iostream>

using std::begin, std::end;

namespace
{
    template<typename R1, typename R2>
    auto comp(const R1& l, const R2& r)
    {
        auto sc = l.size() <=> r.size();
        if (sc != 0) return sc;
        return std::lexicographical_compare_three_way(std::make_reverse_iterator(end(l)), std::make_reverse_iterator(begin(l)), std::make_reverse_iterator(end(r)), std::make_reverse_iterator(begin(r)));
    }

    template<typename R1>
    bool isNull(const R1& r)
    {
        for (auto d : r) if (d) return false;
        return true;
    }

    template<typename R1, typename R2>
    bool add(R1& target, const R2& source)
    {
        auto dest_it = begin(target);
        auto last = end(target);
        bool carry = false;
        for (auto digit : source)
        {
            auto& current = *(dest_it++);
            auto old = current;
            current += digit;
            bool overflow = current < old;
            if (carry) overflow |= !(++current);
            carry = overflow;
        }
        while (carry && dest_it != last)
        {
            auto& current = *(dest_it++);
            carry = !(++current);
        }
        return carry;
    }

    template<typename R1, typename R2>
    bool sub(R1& target, const R2& source)
    {
        auto dest_it = begin(target);
        auto last = end(target);
        bool carry = false;
        for (auto digit : source)
        {
            auto& current = *(dest_it++);
            auto old = current;
            current -= digit;
            bool underflow = current > old;
            if (carry) underflow |= !(current--);
            carry = underflow;
        }
        while (carry && dest_it != last)
        {
            auto& current = *(dest_it++);
            carry = !(current--);
        }
        return carry;
    }

    template<typename T>
    std::pair<T, T> mulDigits(T a, T b)
    {
        const size_t HalfBits = sizeof(T) * CHAR_BIT / 2;
        const T LowMask = (T{ 1 } << HalfBits) - 1;
        std::pair<T, T> result;
        auto& [low, high] = result;
        low = a * b;
        high = (a >> HalfBits) * (b >> HalfBits)
            + (((a & LowMask) * (b & LowMask) >> HalfBits) + (a >> HalfBits) * (b & LowMask) + (a & LowMask) * (b >> HalfBits) >> HalfBits);

        return result;
    }

    template<typename R1, typename R2, typename R3, typename C>
    bool fmc(const R1& left, const R2& right, R3& target, C&& compose)
    {
        auto t_it = begin(target);
        auto t_end = end(target);
        bool overflow = false;
        for (auto l_digit : left)
        {
            auto local_t_it = t_it++;
            for (auto r_digit : right)
            {
                decltype(l_digit) arr[2];
                std::tie(arr[0], arr[1]) = mulDigits(l_digit, r_digit);
                auto subtarget = std::ranges::subrange(local_t_it++, t_end);
                overflow |= compose(subtarget, arr);
            }
        }
        return overflow;
    }

    namespace disambiguity
    {
        template<typename R1, typename R2, typename R3>
        bool fma(const R1& left, const R2& right, R3& target)
        {
            return fmc(left, right, target, [](auto&& l, auto&& r) { return add(std::forward<decltype(l)>(l), std::forward<decltype(r)>(r)); });
        }
    }
    using namespace disambiguity;

    template<typename R1, typename R2, typename R3>
    bool fms(const R1& left, const R2& right, R3& target)
    {
        return fmc(left, right, target, [](auto&& l, auto&& r) { return sub(std::forward<decltype(l)>(l), std::forward<decltype(r)>(r)); });
    }

    template<typename R>
    auto lhs(R& x, size_t s)
    {
        using ValueType = std::remove_cvref_t<decltype(x[0])>;

        if (!s) return ValueType{};
        const size_t bits = sizeof(ValueType) * CHAR_BIT;
        const size_t anti_shift = bits - s;

        ValueType lost = 0;
        for (auto& digit : x)
        {
            ValueType new_lost = digit >> anti_shift;
            digit <<= s;
            digit |= lost;
            lost = new_lost;
        }

        return lost;
    }

    template<typename R>
    auto rhs(R& x, size_t s)
    {
        using ValueType = std::remove_cvref_t<decltype(x[0])>;
        if (!s) return ValueType{};
        const size_t bits = sizeof(ValueType) * CHAR_BIT;
        const size_t anti_shift = bits - s;

        ValueType lost = 0;
        for (auto& digit : std::ranges::reverse_view(x))
        {
            ValueType new_lost = digit << anti_shift;
            digit >>= s;
            digit |= lost;
            lost = new_lost;
        }

        return lost;
    }

    template<typename R>
    auto trim(R& range)
    {
        auto it = std::make_reverse_iterator(end(range));
        auto s = std::make_reverse_iterator(begin(range));
        while (it != s && !(*it)) ++it;
        return std::ranges::subrange(s.base(), it.base());
    }

    template<typename T>
    T wideDivide(T& num_low, T& num_high, T denom)
    {
        const T Factor = static_cast<T>(-1) / denom + ((static_cast<T>(-1) % denom) + 1 == denom);
        T quotient = 0;
        do
        {
            T partial_quotient = num_low / denom + num_high * Factor;
            T num_arr[] = { num_low, num_high };
            fms(std::ranges::single_view(partial_quotient), std::ranges::single_view(denom), num_arr);
            num_low = num_arr[0]; num_high = num_arr[1];
            quotient += partial_quotient;
        } while (num_high || num_low >= denom);
        return quotient;
    }

    template<typename R1, typename T, typename R3>
    void divSimpleDenom(R1& remainder, T denom, R3& quotient)
    {
        using ValueType = std::remove_cvref_t<decltype(remainder[0])>;
        ValueType empty = 0;
        ValueType* previous = &empty;
        auto quo_it = begin(quotient) + remainder.size() - 1;
        for (auto& digit : std::ranges::reverse_view(remainder))
        {
            ValueType partial_quotient = wideDivide(digit, *previous, denom);
            //add(std::ranges::subrange(quo_it, quo_end), std::ranges::single_view(partial_quotient));
            *quo_it = partial_quotient;
            previous = &digit;
            --quo_it;
        }
    }

    template<typename T>
    size_t clz(T x)
    {
        const size_t bits = sizeof(T) * CHAR_BIT;
        T mask = T{ 1 } << (bits - 1);
        size_t result = 0;
        while (mask && !(mask & x))
        {
            ++result;
            mask >>= 1;
        }

        return result;
    }
    
    template<typename T>
    std::pair<T, size_t> maximizing_shift(T val, T secondary={})
    {
        std::pair<T, size_t> result;
        auto& [shifted, shift] = result;
        if(val == static_cast<T>(-1))
        {
            shifted = val >> 1;
            shift = sizeof(T) * CHAR_BIT - 1;
        }
        else
        {
            shift = clz(val);
            T arr[] = {secondary, val};
            lhs(arr, shift);
            shifted = arr[1];
            if(shifted == static_cast<T>(-1))
            {
                shifted >>= 1;
                --shift;
            }
        }
        
        return result;
    }

    namespace disambiguity
    {
        template<typename R1, typename R2, typename R3>
        void div(R1& remainder, const R2& denom, R3& quotient)
        {
            if(denom.size() == 1)
            {
                divSimpleDenom(remainder, denom.front(), quotient);
                return;
            }
            
            using ValueType = std::remove_cvref_t<decltype(remainder[0])>;
            
            auto [simple_denom, shift] = maximizing_shift(denom.back(), denom[denom.size() - 2]);
            bool additional_denom_digit = simple_denom < denom.back();
            ++simple_denom;
            
            auto significant_rem_it = begin(remainder) + denom.size() + (additional_denom_digit ? 0 : -1);
            std::vector<ValueType> rem_copy_buffer(end(remainder) - significant_rem_it + 2);
            std::vector<ValueType> partial_quotient_buffer(rem_copy_buffer.size());
            bool quotient_null;
            do
            {
                std::copy(significant_rem_it - 1, end(remainder), begin(rem_copy_buffer));
                rem_copy_buffer.back() = 0;
                lhs(rem_copy_buffer, shift);
                
                auto rem_copy_untrimmed = std::ranges::subrange(begin(rem_copy_buffer)+1, end(rem_copy_buffer));
                auto rem_copy = trim(rem_copy_untrimmed);
                auto partial_quotient = std::ranges::subrange(begin(partial_quotient_buffer),
                    begin(partial_quotient_buffer) + rem_copy.size());
                
                divSimpleDenom(rem_copy, simple_denom, partial_quotient);

                add(quotient, partial_quotient);
                fms(denom, partial_quotient, remainder);
                
                quotient_null = isNull(partial_quotient);
            } while (!quotient_null);

            while (comp(trim(remainder), trim(denom)) > 0)
            {
                add(quotient, std::ranges::single_view(1));
                sub(remainder, denom);
            }
        }
    }
}

namespace BigNumber
{
    BigUint::BigUint(underlying_t val)
    {
        if(val) _data.push_back(val);
    }
    
    BigUint::underlying_t BigUint::lastDigit() const
    {
        return _data.empty() ? 0 : _data.front();
    }
    
    void BigUint::reset()
    {
        _data.clear();
    }
    
    void BigUint::shrink_to_fit()
    {
        _data.shrink_to_fit();
    }
    
    BigUint fma(const BigUint& l, const BigUint& r, BigUint to)
    {
        to._data.resize(std::max(to._data.size(), l._data.size() + r._data.size()));
        if (disambiguity::fma(l._data, r._data, to._data)) to._data.push_back(1);
        to.shrink();
        return to;
    }
    
    std::pair<BigUint, BigUint> div(BigUint num, const BigUint& denom)
    {
        std::pair<BigUint, BigUint> result;
        auto& [quotient, remainder] = result;
        remainder = std::move(num);
        if(remainder >= denom)
        {
            quotient._data.resize(remainder._data.size() - denom._data.size() + 1);
            disambiguity::div(remainder._data, denom._data, quotient._data);
            remainder.shrink();
            quotient.shrink();
        }
        return result;
    }

    BigUint& BigUint::operator+=(const BigUint& other)
    {
        _data.resize(std::max(_data.size(), other._data.size()));
        if (add(_data, other._data)) _data.push_back(1);

        return *this;
    }

    BigUint& BigUint::operator-=(const BigUint& other)
    {
        _data.resize(std::max(_data.size(), other._data.size()));
        sub(_data, other._data);
        shrink();
        return *this;
    }

    BigUint& BigUint::operator*=(const BigUint& other)
    {
        *this = fma(*this, other, {});

        return *this;
    }

    BigUint& BigUint::operator/=(const BigUint& other)
    {
        *this = div(std::move(*this), other).first;
        
        return *this;
    }
    
    BigUint& BigUint::operator%=(const BigUint& other)
    {
        *this = div(std::move(*this), other).second;
        
        return *this;
    }
    
    void BigUint::shrink()
    {
        _data.resize(trim(_data).size());
    }
    
    BigUint::operator bool() const
    {
        return _data.size();
    }
    
    std::strong_ordering operator<=>(const BigUint& l, const BigUint& r)
    {
        return comp(l._data, r._data);
    }
    
    BigUint operator+(BigUint l, const BigUint& r) { l+=r; return l; }
    BigUint operator-(BigUint l, const BigUint& r) { l-=r; return l; }
    BigUint operator*(BigUint l, const BigUint& r) { l*=r; return l; }
    BigUint operator/(BigUint l, const BigUint& r) { l/=r; return l; }
    BigUint operator%(BigUint l, const BigUint& r) { l%=r; return l; }
    
    std::ostream& operator<<(std::ostream& out, BigUint x)
    {
        if(!x)
        {
            out << 0;
            return out;
        }
        BigUint r;
        std::vector<int> d;
        while(x)
        {
            std::tie(x,r) = div(x, 10);
            d.push_back(r.lastDigit());
        }
        for(int v : std::ranges::reverse_view(d))
        {
            out << v;
        }
        
        return out;
    }
    
    std::istream& operator>>(std::istream& in, BigUint& x)
    {
        using traits = std::istream::traits_type;
        x.reset();
        std::istream::sentry s(in);
        if(!s) return in;
        char c;
        if(!in.get(c) || c < '0' || c > '9')
        {
            in.setstate(std::ios::failbit);
            return in;
        }
        x = c - '0';
        while(in.peek() != traits::eof() && traits::to_char_type(in.peek()) >= '0' && traits::to_char_type(in.peek()) <= '9')
        {
            in.get(c);
            x *= 10;
            x += c - '0';
        }
        
        return in;
    }
}
