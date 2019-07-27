// SPDX-License-Identifier: MIT
// Copyright (c) 2019 Thomas Vanderbruggen <th.vanderbruggen@gmail.com>

#ifndef SCICPP_CORE_UNITS
#define SCICPP_CORE_UNITS

#include "scicpp/core/meta.hpp"

#include <cmath>
#include <cstdio>
#include <numeric>
#include <ratio>
#include <type_traits>

namespace scicpp::units {

// Represents a physical quantity.
// T is the underlying representation type.
// Dim is the dimension of the quantity, that is the rational of primary dimensions.
// Scale is a rational for scaling the quantity.
// Offset mainly for temperature units conversions (ex. K <-> Celsius).

template <typename T, typename Dim, typename Scale, typename Offset>
struct quantity;

// is_quantity

namespace detail {

template <class T>
struct is_quantity : std::false_type {};

template <typename T, typename Dim, typename Scale, typename Offset>
struct is_quantity<quantity<T, Dim, Scale, Offset>> : std::true_type {};

} // namespace detail

template <class T>
inline constexpr bool is_quantity_v = detail::is_quantity<T>::value;

template <class T>
using enable_if_is_quantity = std::enable_if_t<is_quantity_v<T>, int>;

template <class T>
using disable_if_is_quantity = std::enable_if_t<!is_quantity_v<T>, int>;

// common_quantity

template <typename T,
          typename Dim,
          typename Scale1,
          typename Scale2,
          typename Offset1,
          typename Offset2>
struct common_quantity {
    static_assert(std::is_same_v<Offset1, Offset2>,
                  "Do not add/subtract units with different offset, this is "
                  "confusing. Use an explicit call to quantity_cast.");

  private:
    static constexpr auto gcd_num = std::gcd(Scale1::num, Scale2::num);
    static constexpr auto gcd_den = std::gcd(Scale1::den, Scale2::den);
    using CommonScale =
        typename std::ratio<gcd_num,
                            (Scale1::den / gcd_den) * Scale2::den>::type;

  public:
    using type = quantity<T, Dim, CommonScale, Offset1>;
};

template <typename T,
          typename Dim,
          typename Scale1,
          typename Scale2,
          typename Offset1,
          typename Offset2>
using common_quantity_t =
    typename common_quantity<T, Dim, Scale1, Scale2, Offset1, Offset2>::type;

// quantity_cast

template <typename ToQty,
          typename T,
          typename Dim,
          typename Scale,
          typename Offset,
          enable_if_is_quantity<ToQty> = 0>
constexpr auto quantity_cast(const quantity<T, Dim, Scale, Offset> &qty) {
    constexpr auto to_qty_inv_scale =
        T(ToQty::scale::den) / T(ToQty::scale::num);
    constexpr auto qty_scale = T(Scale::num) / T(Scale::den);
    using OffsetDiff = std::ratio_subtract<Offset, typename ToQty::offset>;
    constexpr auto offset_diff = T(OffsetDiff::num) / T(OffsetDiff::den);
    return ToQty(to_qty_inv_scale * (qty_scale * qty.value() + offset_diff));
}

template <typename T,
          typename Dim,
          typename Scale,
          typename Offset = std::ratio<0>>
struct quantity {
    static_assert(meta::is_ratio_v<Dim>);
    static_assert(Dim::num > 0);
    static_assert(Dim::den != 0);
    static_assert(meta::is_ratio_v<Scale>);
    static_assert(Scale::num > 0);
    static_assert(Scale::den != 0);
    static_assert(meta::is_ratio_v<Offset>);
    static_assert(Offset::den != 0);

    using value_type = T;
    using dim = Dim;
    using scale = Scale;
    using offset = Offset;

    // Constructors, destructors, copy

    constexpr quantity() = default;

    quantity(const quantity &) = default;

    explicit constexpr quantity(T value) : m_value(value) {}

    template <typename Scale2, typename Offset2>
    constexpr quantity(const quantity<T, Dim, Scale2, Offset2> &qty)
        : m_value(quantity_cast<quantity>(qty).value()) {}

    ~quantity() = default;

    quantity &operator=(const quantity &) = default;

    // Comparisons

    template <typename RhsScale, typename RhsOffset>
    constexpr auto
    operator==(const quantity<T, Dim, RhsScale, RhsOffset> &rhs) const {
        constexpr auto lhs_factor =
            T(Scale::num * RhsScale::den * Offset::den * RhsOffset::den);
        constexpr auto lhs_offset =
            T(Offset::num * RhsScale::den * Scale::den * RhsOffset::den);
        constexpr auto rhs_factor =
            T(RhsScale::num * Scale::den * RhsOffset::den * Offset::den);
        constexpr auto rhs_offset =
            T(RhsOffset::num * Scale::den * RhsScale::den * Offset::den);

        return lhs_factor * m_value + lhs_offset ==
               rhs_factor * rhs.value() + rhs_offset;
    }

    template <typename RhsScale, typename RhsOffset>
    constexpr auto
    operator!=(const quantity<T, Dim, RhsScale, RhsOffset> &rhs) const {
        return !(*this == rhs);
    }

    template <typename RhsScale, typename RhsOffset>
    constexpr auto
    operator<(const quantity<T, Dim, RhsScale, RhsOffset> &rhs) const {
        constexpr auto lhs_factor =
            T(Scale::num * RhsScale::den * Offset::den * RhsOffset::den);
        constexpr auto lhs_offset =
            T(Offset::num * RhsScale::den * Scale::den * RhsOffset::den);
        constexpr auto rhs_factor =
            T(RhsScale::num * Scale::den * RhsOffset::den * Offset::den);
        constexpr auto rhs_offset =
            T(RhsOffset::num * Scale::den * RhsScale::den * Offset::den);

        return lhs_factor * m_value + lhs_offset <
               rhs_factor * rhs.value() + rhs_offset;
    }

    template <typename RhsScale, typename RhsOffset>
    constexpr auto
    operator<=(const quantity<T, Dim, RhsScale, RhsOffset> &rhs) const {
        return !(rhs < *this);
    }

    template <typename RhsScale, typename RhsOffset>
    constexpr auto
    operator>(const quantity<T, Dim, RhsScale, RhsOffset> &rhs) const {
        return rhs < *this;
    }

    template <typename RhsScale, typename RhsOffset>
    constexpr auto
    operator>=(const quantity<T, Dim, RhsScale, RhsOffset> &rhs) const {
        return !(*this < rhs);
    }

    // Arithmetic

    constexpr auto &operator+=(const quantity &rhs) {
        m_value += rhs.m_value;
        return *this;
    }

    constexpr auto &operator-=(const quantity &rhs) {
        m_value -= rhs.m_value;
        return *this;
    }

    constexpr auto operator+() const { return *this; }
    constexpr auto operator-() const { return quantity(-m_value); }

    constexpr const quantity operator++(int) { return quantity(m_value++); }

    constexpr quantity &operator++() {
        ++m_value;
        return *this;
    }

    constexpr const quantity operator--(int) { return quantity(m_value--); }

    constexpr quantity &operator--() {
        --m_value;
        return *this;
    }

    template <typename RhsScale, typename RhsOffset>
    constexpr auto
    operator+(const quantity<T, Dim, RhsScale, RhsOffset> &rhs) const {
        using cq =
            common_quantity_t<T, Dim, Scale, RhsScale, Offset, RhsOffset>;
        return cq(cq(*this).value() + cq(rhs).value());
    }

    template <typename RhsScale, typename RhsOffset>
    constexpr auto
    operator-(const quantity<T, Dim, RhsScale, RhsOffset> &rhs) const {
        using cq =
            common_quantity_t<T, Dim, Scale, RhsScale, Offset, RhsOffset>;
        return cq(cq(*this).value() - cq(rhs).value());
    }

    template <typename RhsDim, typename RhsScale>
    constexpr auto operator*(const quantity<T, RhsDim, RhsScale> &rhs) const {
        return quantity<T,
                        std::ratio_multiply<Dim, RhsDim>,
                        std::ratio_multiply<Scale, RhsScale>>(m_value *
                                                              rhs.value());
    }

    template <typename RhsDim, typename RhsScale>
    constexpr auto operator/(const quantity<T, RhsDim, RhsScale> &rhs) const {
        return quantity<T,
                        std::ratio_divide<Dim, RhsDim>,
                        std::ratio_divide<Scale, RhsScale>>(m_value /
                                                            rhs.value());
    }

    constexpr auto inv() const {
        return quantity<T,
                        std::ratio_divide<std::ratio<1>, Dim>,
                        std::ratio_divide<std::ratio<1>, Scale>>(T{1} /
                                                                 m_value);
    }

    constexpr auto value() const { return m_value; }

  private:
    T m_value;
};

template <typename T1,
          typename T2,
          typename Dim,
          typename Scale,
          typename Offset,
          meta::disable_if_iterable<T1> = 0>
constexpr auto operator*(T1 factor,
                         const quantity<T2, Dim, Scale, Offset> &rhs) {
    using T = decltype(std::declval<T1>() * std::declval<T2>());
    return quantity<T, Dim, Scale, Offset>(T(factor) * T(rhs.value()));
}

template <typename T1,
          typename T2,
          typename Dim,
          typename Scale,
          typename Offset,
          meta::disable_if_iterable<T2> = 0>
constexpr auto operator*(const quantity<T1, Dim, Scale, Offset> &rhs,
                         T2 factor) {
    return factor * rhs;
}

template <typename T1,
          typename T2,
          typename Dim,
          typename Scale,
          typename Offset,
          meta::disable_if_iterable<T1> = 0>
constexpr auto operator/(T1 factor,
                         const quantity<T2, Dim, Scale, Offset> &rhs) {
    return factor * rhs.inv();
}

template <typename T1,
          typename T2,
          typename Dim,
          typename Scale,
          typename Offset,
          meta::disable_if_iterable<T2> = 0>
constexpr auto operator/(const quantity<T1, Dim, Scale, Offset> &rhs,
                         T2 factor) {
    return rhs * (T2{1} / factor);
}

template <typename T, typename Dim, typename Scale>
bool isnan(const quantity<T, Dim, Scale> &q) {
    return std::isnan(q.value());
}

template <typename T>
auto value(T x) {
    if constexpr (is_quantity_v<T>) {
        return x.value();
    } else {
        return x;
    }
}

template <typename T, typename Dim, typename Scale, typename Offset>
void print(const quantity<T, Dim, Scale, Offset> &q) {
    std::printf("%f x (%li / %li) + (%li / %li) = %f [%li, %li]\n",
                q.value(),
                Scale::num,
                Scale::den,
                Offset::num,
                Offset::den,
                q.value() * (T(Scale::num) / T(Scale::den)) +
                    T(Offset::num) / T(Offset::den),
                Dim::num,
                Dim::den);
}

// representation_t: Get the representation type of a quantity

namespace detail {

template <class T>
struct representation_type_impl {
    using type = T;
};

template <typename T, typename Dim, typename Scale, typename Offset>
struct representation_type_impl<quantity<T, Dim, Scale, Offset>> {
    using type = T;
};

} // namespace detail

template <class T>
using representation_t = typename detail::representation_type_impl<T>::type;

// Quantities type operations

template <typename Quantity1, typename Quantity2>
using quantity_multiply =
    decltype(std::declval<Quantity1>() * std::declval<Quantity2>());

template <typename Quantity1, typename Quantity2>
using quantity_divide =
    decltype(std::declval<Quantity1>() / std::declval<Quantity2>());

template <typename Quantity>
using quantity_invert = decltype(std::declval<Quantity>().inv());

namespace primary_flags {

// Each base quantity is identified by a prime number.
// using Dimensionless = std::ratio<1>;

// SI
using Length = std::ratio<2>;
using Time = std::ratio<3>;
using Mass = std::ratio<5>;
using ElectricCurrent = std::ratio<7>;
using Temperature = std::ratio<11>;
using AmountOfSubstance = std::ratio<13>;
using LuminousIntensity = std::ratio<17>;

// Data
using DataQuantity = std::ratio<19>;

} // namespace primary_flags

#define SCICPP_CORE_UNITS_SET_LITERAL(quantity, literal, scale)                \
    constexpr auto operator""##literal(long double x) {                        \
        return quantity<double, scale>{static_cast<double>(x)};                \
    }                                                                          \
                                                                               \
    constexpr auto operator""##literal(unsigned long long x) {                 \
        return quantity<double, scale>{static_cast<double>(x)};                \
    }

#define SCICPP_CORE_UNITS_SET_LITERAL_RATIO(quantity, literal, num, den)       \
    constexpr auto operator""##literal(long double x) {                        \
        return quantity<double, std::ratio<(num), (den)>>{                     \
            static_cast<double>(x)};                                           \
    }                                                                          \
                                                                               \
    constexpr auto operator""##literal(unsigned long long x) {                 \
        return quantity<double, std::ratio<(num), (den)>>{                     \
            static_cast<double>(x)};                                           \
    }

#define SCICPP_CORE_UNITS_SET_LITERAL_RATIO2(                                  \
    quantity, literal, num1, den1, num2, den2)                                 \
    constexpr auto operator""##literal(long double x) {                        \
        return quantity<double,                                                \
                        std::ratio<(num1), (den1)>,                            \
                        std::ratio<(num2), (den2)>>{static_cast<double>(x)};   \
    }                                                                          \
                                                                               \
    constexpr auto operator""##literal(unsigned long long x) {                 \
        return quantity<double,                                                \
                        std::ratio<(num1), (den1)>,                            \
                        std::ratio<(num2), (den2)>>{static_cast<double>(x)};   \
    }

#define SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(quantity, base_unit)           \
    template <typename T = double>                                             \
    using atto##base_unit = quantity<T, std::atto>;                            \
    template <typename T = double>                                             \
    using femto##base_unit = quantity<T, std::femto>;                          \
    template <typename T = double>                                             \
    using pico##base_unit = quantity<T, std::pico>;                            \
    template <typename T = double>                                             \
    using nano##base_unit = quantity<T, std::nano>;                            \
    template <typename T = double>                                             \
    using micro##base_unit = quantity<T, std::micro>;                          \
    template <typename T = double>                                             \
    using milli##base_unit = quantity<T, std::milli>;                          \
    template <typename T = double>                                             \
    using centi##base_unit = quantity<T, std::centi>;                          \
    template <typename T = double>                                             \
    using deci##base_unit = quantity<T, std::deci>;                            \
    template <typename T = double>                                             \
    using base_unit = quantity<T>;                                             \
    template <typename T = double>                                             \
    using deca##base_unit = quantity<T, std::deca>;                            \
    template <typename T = double>                                             \
    using hecto##base_unit = quantity<T, std::hecto>;                          \
    template <typename T = double>                                             \
    using kilo##base_unit = quantity<T, std::kilo>;                            \
    template <typename T = double>                                             \
    using mega##base_unit = quantity<T, std::mega>;                            \
    template <typename T = double>                                             \
    using giga##base_unit = quantity<T, std::giga>;                            \
    template <typename T = double>                                             \
    using tera##base_unit = quantity<T, std::tera>;                            \
    template <typename T = double>                                             \
    using peta##base_unit = quantity<T, std::peta>;                            \
    template <typename T = double>                                             \
    using exa##base_unit = quantity<T, std::exa>;

// Primary quantities

// ----------------------------------------------------------------------------
// Length
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using length = quantity<T, primary_flags::Length, Scale>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(length, meter)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(length, _fm, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(length, _pm, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(length, _nm, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(length, _um, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(length, _mm, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(length, _cm, std::centi)
SCICPP_CORE_UNITS_SET_LITERAL(length, _m, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(length, _km, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(length, _Mm, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(length, _Gm, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(length, _Tm, std::tera)

SCICPP_CORE_UNITS_SET_LITERAL_RATIO(length, _angstrom, 1, 10000000000)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(length, _in, 254, 10000)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(length, _ft, 3048, 10000)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(length, _yd, 9144, 10000)

} // namespace literals

// ----------------------------------------------------------------------------
// Time
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using time = quantity<T, primary_flags::Time, Scale>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(time, second)

template <typename T = double>
using minute = time<T, std::ratio<60>>;

template <typename T = double>
using hour = time<T, std::ratio<3600>>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(time, _fs, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(time, _ps, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(time, _ns, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(time, _us, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(time, _ms, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(time, _s, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(time, _min, 60, 1)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(time, _h, 3600, 1)

} // namespace literals

// ----------------------------------------------------------------------------
// Mass
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using mass = quantity<T, primary_flags::Mass, Scale>;

template <typename T = double>
using gram = mass<T, std::milli>;

template <typename T = double>
using kilogram = mass<T>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(mass, _ng, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _ug, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _mg, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _g, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _kg, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _ton, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _kton, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _Mton, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _Gton, std::tera)
SCICPP_CORE_UNITS_SET_LITERAL(mass, _Tton, std::exa)

} // namespace literals

// ----------------------------------------------------------------------------
// Electric current
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using electric_current = quantity<T, primary_flags::ElectricCurrent, Scale>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(electric_current, ampere)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _fA, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _pA, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _nA, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _uA, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _mA, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _A, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _kA, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _MA, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _GA, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(electric_current, _TA, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Temperature
// ----------------------------------------------------------------------------

template <typename T,
          typename Scale = std::ratio<1>,
          typename Offset = std::ratio<0>>
using temperature = quantity<T, primary_flags::Temperature, Scale, Offset>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(temperature, kelvin)

template <typename T = double>
using celsius = temperature<T, std::ratio<1>, std::ratio<27315, 100>>;

template <typename T = double>
using fahrhenheit = temperature<T, std::ratio<5, 9>, std::ratio<45967, 180>>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(temperature, _fK, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(temperature, _pK, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(temperature, _nK, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(temperature, _uK, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(temperature, _mK, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(temperature, _K, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO2(temperature, _degC, 1, 1, 27315, 100)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO2(temperature, _degF, 5, 9, 45967, 180)

} // namespace literals

// ----------------------------------------------------------------------------
// Amount of substance
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using amount_of_substance =
    quantity<T, primary_flags::AmountOfSubstance, Scale>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(amount_of_substance, mole)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(amount_of_substance, _fmol, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(amount_of_substance, _pmol, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(amount_of_substance, _nmol, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(amount_of_substance, _umol, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(amount_of_substance, _mmol, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(amount_of_substance, _mol, std::ratio<1>)

} // namespace literals

// ----------------------------------------------------------------------------
// Luminous intensity
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using luminous_intensity = quantity<T, primary_flags::LuminousIntensity, Scale>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(luminous_intensity, candela)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _fCd, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _pCd, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _nCd, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _uCd, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _mCd, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _Cd, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _kCd, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _MCd, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(luminous_intensity, _GCd, std::giga)

} // namespace literals

// ----------------------------------------------------------------------------
// Data quantity
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using data_quantity = quantity<T, primary_flags::DataQuantity, Scale>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(data_quantity, _b, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(data_quantity, _B, 8, 1)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(data_quantity, _kiB, 8LL * 1024LL, 1)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(data_quantity,
                                    _MiB,
                                    8LL * 1024LL * 1024LL,
                                    1)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(data_quantity,
                                    _GiB,
                                    8LL * 1024LL * 1024LL * 1024LL,
                                    1)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(data_quantity,
                                    _TiB,
                                    8LL * 1024LL * 1024LL * 1024LL * 1024LL,
                                    1)

} // namespace literals

// Derived quantities

// ----------------------------------------------------------------------------
// Speed
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using speed = quantity_divide<length<T, Scale>, time<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(speed, meter_per_second)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(speed, _fm_per_s, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(speed, _pm_per_s, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(speed, _nm_per_s, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(speed, _um_per_s, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(speed, _mm_per_s, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(speed, _m_per_s, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(speed, _km_per_s, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(speed, _km_per_h, 3600, 1000)

} // namespace literals

// ----------------------------------------------------------------------------
// Acceleration
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using acceleration = quantity_divide<speed<T, Scale>, time<T>>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(acceleration, _fm_per_s2, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(acceleration, _pm_per_s2, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(acceleration, _nm_per_s2, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(acceleration, _um_per_s2, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(acceleration, _mm_per_s2, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(acceleration, _m_per_s2, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(acceleration, _km_per_s2, std::kilo)

} // namespace literals

// ----------------------------------------------------------------------------
// area
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using area = quantity_multiply<length<T, Scale>, length<T>>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(area, _m2, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(area, _km2, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(area, _ha, 10000, 1)

} // namespace literals

// ----------------------------------------------------------------------------
// Volume
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using volume = quantity_multiply<area<T, Scale>, length<T>>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(volume, _m3, std::ratio<1>)

} // namespace literals

// ----------------------------------------------------------------------------
// Force
// ----------------------------------------------------------------------------

// Force = Mass x Acceleration
template <typename T, typename Scale = std::ratio<1>>
using force = quantity_multiply<mass<T, Scale>, acceleration<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(force, newton)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(force, _fN, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(force, _pN, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(force, _nN, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(force, _uN, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(force, _mN, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(force, _N, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(force, _kN, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(force, _MN, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(force, _GN, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(force, _TN, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Power
// ----------------------------------------------------------------------------

// Power = Force x Speed
template <typename T, typename Scale = std::ratio<1>>
using power = quantity_multiply<force<T, Scale>, speed<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(power, watt)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(power, _fW, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(power, _pW, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(power, _nW, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(power, _uW, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(power, _mW, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(power, _W, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(power, _kW, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(power, _MW, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(power, _GW, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(power, _TW, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Energy
// ----------------------------------------------------------------------------

// Energy = Power x Time
template <typename T, typename Scale = std::ratio<1>>
using energy = quantity_multiply<power<T, Scale>, time<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(energy, joule)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(energy, _fJ, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _pJ, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _nJ, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _uJ, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _mJ, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _J, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _kJ, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _MJ, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _GJ, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(energy, _TJ, std::tera)

SCICPP_CORE_UNITS_SET_LITERAL_RATIO(energy, _kWh, 3600000, 1)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(energy, _cal, 41855, 10000)

} // namespace literals

// ----------------------------------------------------------------------------
// Pressure
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using pressure = quantity_divide<force<T, Scale>, area<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(pressure, pascal)

template <typename T = double>
using bar = pressure<T, std::ratio<100000>>;

template <typename T = double>
using mmHg = pressure<T, std::ratio<101325, 760>>;

template <typename T = double>
using torr = mmHg<T>;

template <typename T = double>
using psi = pressure<T, std::ratio<689476, 100>>;

template <typename T = double>
using atm = pressure<T, std::ratio<101325>>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(pressure, _fPa, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _pPa, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _nPa, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _uPa, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _mPa, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _Pa, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _kPa, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _MPa, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _GPa, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(pressure, _TPa, std::tera)

SCICPP_CORE_UNITS_SET_LITERAL_RATIO(pressure, _bar, 100000, 1)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(pressure, _mmHg, 101325, 760)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(pressure, _torr, 101325, 760)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(pressure, _psi, 689476, 100)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(pressure, _atm, 101325, 1)

} // namespace literals

// ----------------------------------------------------------------------------
// Frequency
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using frequency = quantity_invert<time<T, Scale>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(frequency, hertz)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(frequency, _fHz, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _pHz, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _nHz, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _uHz, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _mHz, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _Hz, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _kHz, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _MHz, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _GHz, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(frequency, _THz, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Electrical resistance
// ----------------------------------------------------------------------------

// Resistance = Power / Current^2
template <typename T, typename Scale = std::ratio<1>>
using electrical_resistance = quantity_divide<
    power<T, Scale>,
    quantity_multiply<electric_current<T>, electric_current<T>>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(electrical_resistance, ohm)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _fOhm, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _pOhm, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _nOhm, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _uOhm, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _mOhm, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _Ohm, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _kOhm, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _MOhm, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _GOhm, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_resistance, _TOhm, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Electrical conductivity
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using electrical_conductivity =
    quantity_invert<electrical_resistance<T, Scale>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(electrical_conductivity, siemens)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _fS, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _pS, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _nS, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _uS, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _mS, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _S, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _kS, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _MS, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _GS, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(electrical_conductivity, _TS, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Electric potential
// ----------------------------------------------------------------------------

// Voltage = Power / Current
template <typename T, typename Scale = std::ratio<1>>
using electric_potential =
    quantity_divide<power<T, Scale>, electric_current<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(electric_potential, volt)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _fV, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _pV, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _nV, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _uV, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _mV, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _V, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _kV, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _MV, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _GV, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(electric_potential, _TV, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Capacitance
// ----------------------------------------------------------------------------

// Capacitance = Time / Resistance
template <typename T, typename Scale = std::ratio<1>>
using capacitance = quantity_divide<time<T, Scale>, electrical_resistance<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(capacitance, farad)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _fF, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _pF, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _nF, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _uF, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _mF, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _F, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _kF, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _MF, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _GF, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(capacitance, _TF, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Electric charge
// ----------------------------------------------------------------------------

// Charge = Capacitance x Voltage
template <typename T, typename Scale = std::ratio<1>>
using electric_charge =
    quantity_multiply<capacitance<T, Scale>, electric_potential<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(electric_charge, coulomb)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _fC, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _pC, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _nC, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _uC, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _mC, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _C, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _kC, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _MC, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _GC, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(electric_charge, _TC, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Inductance
// ----------------------------------------------------------------------------

// Inductance = Voltage x Time / Current
template <typename T, typename Scale = std::ratio<1>>
using inductance =
    quantity_divide<quantity_multiply<electric_potential<T, Scale>, time<T>>,
                    electric_current<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(inductance, henry)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(inductance, _fH, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _pH, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _nH, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _uH, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _mH, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _H, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _kH, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _MH, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _GH, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(inductance, _TH, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Magnetic flux density
// ----------------------------------------------------------------------------

// Magnetic field = Force / (Charge x Speed)
template <typename T, typename Scale = std::ratio<1>>
using magnetic_flux_density =
    quantity_divide<force<T, Scale>,
                    quantity_multiply<electric_charge<T>, speed<T>>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(magnetic_flux_density, tesla)

template <typename T = double>
using gauss = magnetic_flux_density<T, std::ratio<1, 10000>>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux_density, _fT, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux_density, _pT, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux_density, _nT, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux_density, _uT, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux_density, _mT, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux_density, _T, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL_RATIO(magnetic_flux_density, _G, 1, 10000)

} // namespace literals

// ----------------------------------------------------------------------------
// Magnetic flux
// ----------------------------------------------------------------------------

// Magnetic field = Force / (Charge x Speed)
template <typename T, typename Scale = std::ratio<1>>
using magnetic_flux =
    quantity_multiply<magnetic_flux_density<T, Scale>, area<T>>;

SCICPP_CORE_UNITS_DEFINE_PREFIXES_ALIAS(magnetic_flux, weber)

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _fWb, std::femto)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _pWb, std::pico)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _nWb, std::nano)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _uWb, std::micro)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _mWb, std::milli)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _Wb, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _kWb, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _MWb, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _GWb, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(magnetic_flux, _TWb, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Data rate
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using data_rate = quantity_divide<data_quantity<T, Scale>, time<T>>;

namespace literals {

SCICPP_CORE_UNITS_SET_LITERAL(data_rate, _bps, std::ratio<1>)
SCICPP_CORE_UNITS_SET_LITERAL(data_rate, _kbps, std::kilo)
SCICPP_CORE_UNITS_SET_LITERAL(data_rate, _Mbps, std::mega)
SCICPP_CORE_UNITS_SET_LITERAL(data_rate, _Gbps, std::giga)
SCICPP_CORE_UNITS_SET_LITERAL(data_rate, _Tbps, std::tera)

} // namespace literals

// ----------------------------------------------------------------------------
// Chemical potential
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using chemical_potential =
    quantity_divide<energy<T, Scale>, amount_of_substance<T>>;

// ----------------------------------------------------------------------------
// Thermal resistance
// ----------------------------------------------------------------------------

template <typename T, typename Scale = std::ratio<1>>
using thermal_resistance = quantity_divide<temperature<T, Scale>, power<T>>;

} // namespace scicpp::units

#endif // SCICPP_CORE_UNITS