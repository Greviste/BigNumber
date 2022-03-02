#ifndef BIGNUMBER_BIGNUMBER_H
#define BIGNUMBER_BIGNUMBER_H

#include <vector>
#include <utility>
#include <compare>
#include <cstdint>
#include <ranges>
#include <initializer_list>
#include <iostream>

namespace BigNumber
{
    struct BigUint
    {
    public:
        using underlying_t = std::uintmax_t;
        
        BigUint(underlying_t val=0);
        
        BigUint(std::initializer_list<underlying_t> l) :BigUint(std::ranges::views::all(l)) {}
        
        template<std::ranges::range R>
        explicit BigUint(const R& r)
            :_data(begin(r), end(r))
        {
            shrink();
        }
        
        underlying_t lastDigit() const;
        
        //Return type models borrowed_range, sized_range, common_range, bidirectional_range
        //The first element is the least significant digit
        auto& digits() const
        {
            return _data;
        }
        
        void reset();
        
        friend BigUint fma(const BigUint&, const BigUint&, BigUint);
        friend std::pair<BigUint, BigUint> div(BigUint, const BigUint&); //Returns quotient and remainder
        BigUint& operator+=(const BigUint&);
        BigUint& operator-=(const BigUint&); //Value is unspecified in case of underflow
        BigUint& operator*=(const BigUint&);
        BigUint& operator/=(const BigUint&);
        BigUint& operator%=(const BigUint&);
        explicit operator bool() const;
        friend std::strong_ordering operator<=>(const BigUint&, const BigUint&);
    private:
        void shrink();
        
        std::vector<underlying_t> _data;
    };
    
    BigUint operator+(BigUint, const BigUint&);
    BigUint operator-(BigUint, const BigUint&);
    BigUint operator*(BigUint, const BigUint&);
    BigUint operator/(BigUint, const BigUint&);
    BigUint operator%(BigUint, const BigUint&);
    std::ostream& operator<<(std::ostream&, BigUint);
    std::istream& operator>>(std::istream&, BigUint&);
}

#endif
