// SPDX-License-Identifier: MIT
// Copyright (c) 2019 Thomas Vanderbruggen <th.vanderbruggen@gmail.com>

#ifndef SCICPP_CORE_ARITHMETIC
#define SCICPP_CORE_ARITHMETIC

#include "scicpp/core/meta.hpp"

#include <array>
#include <cstdint>
#include <numeric>
#include <ratio>
#include <utility>

namespace scicpp::arithmetic {

//---------------------------------------------------------------------------------
// Compile-time log2
// https://hbfs.wordpress.com/2016/03/22/log2-with-c-metaprogramming/
//---------------------------------------------------------------------------------

constexpr intmax_t ct_log2(intmax_t num) {
    return (num < 2) ? 1 : 1 + ct_log2(num >> 1);
}

//---------------------------------------------------------------------------------
// Compile-time power
//https://stackoverflow.com/questions/27270541/is-there-no-built-in-way-to-compute-a-power-at-compile-time-in-c
//---------------------------------------------------------------------------------

template <typename T>
constexpr T power(T a, intmax_t n) {
    if (n == 0) {
        return 1;
    }

    const auto p = power(a, n / 2);
    return p * p * (n % 2 == 0 ? 1 : a);
}

//---------------------------------------------------------------------------------
// Compile-time prime factor decomposition
//---------------------------------------------------------------------------------

template <intmax_t n>
struct prime_factors {
    using prime_t = intmax_t;
    using pow_t = intmax_t;

    // https://stackoverflow.com/questions/11924249/finding-prime-factors
    constexpr prime_factors() {
        intmax_t z = 2;
        intmax_t num = n;

        while (z * z <= num) {
            if (num % z == 0) {
                append(z);
                num /= z;
            } else {
                z++;
            }
        }

        if (num > 1) {
            append(num);
        }
    }

    constexpr auto values() const { return factors; }

  private:
    constexpr void append(prime_t prime) {
        std::size_t i = 0;

        for (; i < factors.size(); ++i) {
            if (factors[i].first == 0) {
                break;
            }

            if (factors[i].first == prime) {
                ++factors[i].second;
                return;
            }
        }

        // prime not in array
        factors[i].first = prime;
        factors[i].second = 1;
    }

    // There is at most log2(n) prime factors
    std::array<std::pair<prime_t, pow_t>, std::size_t(ct_log2(n))> factors{};
};

//---------------------------------------------------------------------------------
// Root of a ratio
//---------------------------------------------------------------------------------

// Root of a rational number (Q/P)^(1/R)
// where
// - Q/P is an irreducible fraction
// - Q and P are not perfect Rth power

template <typename Ratio, intmax_t Root = 1>
struct root_ratio {
    static_assert(meta::is_ratio_v<Ratio>);
    static_assert(Ratio::num > 0);
    static_assert(Ratio::den > 0);
    static_assert(Root > 0);

    using ratio = Ratio;

    static constexpr auto num = Ratio::num;
    static constexpr auto den = Ratio::den;
    static constexpr auto root = Root;
};

namespace detail {

template <typename PrimeFactors>
constexpr bool is_exact_root(PrimeFactors factors, intmax_t root) {
    for (const auto &factor : factors) {
        if (factor.second % root != 0) {
            return false;
        }
    }

    return true;
}

template <typename PrimeFactors>
constexpr auto compute_root(PrimeFactors factors, intmax_t root) {
    intmax_t res = 1;

    for (const auto &factor : factors) {
        res *= power(factor.first, factor.second / root);
    }

    return res;
}

template <typename RR, intmax_t Root = 1>
constexpr auto root_ratio_root_impl() {
    static_assert(Root > 0);

    constexpr auto Q = RR::num;
    constexpr auto P = RR::den;
    constexpr auto R = RR::root * Root;

    // If Q and P are exact roots of R (Q=q^R and P=p^R),
    // we can reduce the dimension to a ratio only:
    // (Q/P)^(1/R) = (q^R / p^R)^(1/R) = q/p
    constexpr auto q_factors = prime_factors<Q>().values();
    constexpr auto p_factors = prime_factors<P>().values();

    if constexpr (is_exact_root(q_factors, R) && is_exact_root(p_factors, R)) {
        return root_ratio<std::ratio<compute_root(q_factors, R),
                                     compute_root(p_factors, R)>>{};
    } else {
        return root_ratio<std::ratio<Q, P>, R>{};
    }
}

template <typename RR1, typename RR2>
constexpr auto root_ratio_multiply_impl() {
    constexpr auto R1 = RR1::root;
    constexpr auto R2 = RR2::root;

    if constexpr (R1 == 1 && R2 == 1) {
        return root_ratio<
            std::ratio_multiply<typename RR1::ratio, typename RR2::ratio>>{};
    } else {
        constexpr auto Q1 = RR1::num;
        constexpr auto Q2 = RR2::num;
        constexpr auto P1 = RR1::den;
        constexpr auto P2 = RR2::den;

        // Reduce under common root:
        // (Q1/P1)^(1/R1) x (Q2/P2)^(1/R2) = (Q/P)^(1/R)
        // where:
        // Q = Q1^r2 x Q2^r1
        // P = P1^r2 x P2^r1
        // R = G r1 r2
        // and
        // G = gcd(R1, R2),
        // Ri = G ri

        constexpr auto G = std::gcd(R1, R2);
        constexpr auto r1 = R1 / G;
        constexpr auto r2 = R2 / G;

        constexpr auto Q = power(Q1, r2) * power(Q2, r1);
        constexpr auto P = power(P1, r2) * power(P2, r1);
        constexpr auto R = G * r1 * r2;

        // Compute the root to reduce the root_ratio if exact root
        return root_ratio_root_impl<root_ratio<std::ratio<Q, P>, R>>();
    }
}

template <typename RR1, typename RR2>
constexpr auto root_ratio_divide_impl() {
    using RR2Inv = root_ratio<std::ratio<RR2::den, RR2::num>, RR2::root>;
    return root_ratio_multiply_impl<RR1, RR2Inv>();
}

} // namespace detail

template <typename RR, intmax_t Root = 1>
using root_ratio_root = decltype(detail::root_ratio_root_impl<RR, Root>());

template <typename RR1, typename RR2>
using root_ratio_multiply =
    decltype(detail::root_ratio_multiply_impl<RR1, RR2>());

template <typename RR1, typename RR2>
using root_ratio_divide = decltype(detail::root_ratio_divide_impl<RR1, RR2>());

template <typename RR1, typename RR2>
struct common_root_ratio {
    static constexpr auto Q1_R2 = power(RR1::num, RR2::root);
    static constexpr auto Q2_R1 = power(RR2::num, RR1::root);
    static constexpr auto P1_R2 = power(RR1::den, RR2::root);
    static constexpr auto P2_R1 = power(RR2::den, RR1::root);
    static constexpr auto R = RR1::root * RR2::root;

    static constexpr auto G_Q = std::gcd(Q1_R2, Q2_R1);
    static constexpr auto G_P = std::gcd(Q2_R1, P2_R1);

    static constexpr auto Q = G_Q;
    static constexpr auto P = (P1_R2 / G_P) * P2_R1;

    using type = root_ratio_root<root_ratio<std::ratio<Q, P>, R>>;
};

template <typename RR1, typename RR2>
using common_root_ratio_t = typename common_root_ratio<RR1, RR2>::type;

} // namespace scicpp::arithmetic

#endif // SCICPP_CORE_ARITHMETIC