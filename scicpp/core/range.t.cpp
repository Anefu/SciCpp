// SPDX-License-Identifier: MIT
// Copyright (c) 2019 Thomas Vanderbruggen <th.vanderbruggen@gmail.com>

#include "range.hpp"

#include "scicpp/core/equal.hpp"
#include "scicpp/core/numeric.hpp"
#include "scicpp/core/print.hpp"
#include "scicpp/core/units.hpp"

namespace scicpp {

TEST_CASE("linspace") {
    SECTION("std::array") {
        REQUIRE(linspace<0>(2., 3.).empty());
        REQUIRE(almost_equal(linspace<1>(2., 3.), {2.}));
        REQUIRE(almost_equal(linspace<5>(2., 3.), {2., 2.25, 2.5, 2.75, 3.}));
        REQUIRE(almost_equal(linspace<5>(3., 2.), {3., 2.75, 2.5, 2.25, 2.}));
    }

    SECTION("std::array with quantity") {
        using namespace units::literals;
        REQUIRE(linspace<0>(2._m, 3._m).empty());
        REQUIRE(almost_equal(linspace<1>(2._m, 3._m), {2._m}));
        REQUIRE(almost_equal(linspace<5>(2._m, 3._m),
                             {2._m, 2.25_m, 2.5_m, 2.75_m, 3._m}));
        REQUIRE(almost_equal(linspace<5>(3._kg, 2._kg),
                             {3._kg, 2.75_kg, 2.5_kg, 2.25_kg, 2._kg}));
    }

    SECTION("std::vector") {
        REQUIRE(linspace(2., 3., 0).empty());
        REQUIRE(almost_equal(linspace(2., 3., 1), {2.}));
        REQUIRE(almost_equal(linspace(2., 3., 5), {2., 2.25, 2.5, 2.75, 3.}));
        REQUIRE(almost_equal(linspace(3., 2., 5), {3., 2.75, 2.5, 2.25, 2.}));
    }

    SECTION("std::vector with quantity") {
        using namespace units::literals;
        REQUIRE(linspace(2._m, 3._m, 0).empty());
        REQUIRE(almost_equal(linspace(2._m, 3._m, 1), {2._m}));
        REQUIRE(almost_equal(linspace(2._m, 3._m, 5),
                             {2._m, 2.25_m, 2.5_m, 2.75_m, 3._m}));
        REQUIRE(almost_equal(linspace(3._kg, 2._kg, 5),
                             {3._kg, 2.75_kg, 2.5_kg, 2.25_kg, 2._kg}));
    }
}

TEST_CASE("logspace std::array") {
    SECTION("std::array") {
        REQUIRE(logspace<0>(2., 3.).empty());
        REQUIRE(almost_equal(logspace<1>(2., 3.), {100.}));
        REQUIRE(almost_equal(
            logspace<4>(2., 3.),
            {100., 215.44346900318845, 464.15888336127773, 1000.}));
    }

    SECTION("std::vector") {
        REQUIRE(logspace(2., 3., 0).empty());
        REQUIRE(almost_equal(logspace(2., 3., 1), {100.}));
        REQUIRE(almost_equal(
            logspace(2., 3., 4),
            {100., 215.44346900318845, 464.15888336127773, 1000.}));
    }
}

TEST_CASE("arange std::vector") {
    REQUIRE(arange(2., 3., 10.).empty());
    REQUIRE(almost_equal(arange(2., 3., 0.75), {2.}));
    REQUIRE(almost_equal(arange(0., 3.), {0., 1., 2.}));
    REQUIRE(arange(1., 0., 0.5).empty());
    REQUIRE(arange(0., 3., -1.).empty());
    REQUIRE(almost_equal(arange(3., 0., -1.), {3., 2., 1.}));
}

} // namespace scicpp
