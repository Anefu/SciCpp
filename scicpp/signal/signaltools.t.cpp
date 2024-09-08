#include "signaltools.hpp"

#include "scicpp/core/equal.hpp"
#include "scicpp/signal/filter_design.hpp"
#include "scicpp/core/numeric.hpp"
#include "scicpp/core/print.hpp"

namespace scicpp {
    namespace signal {
        TEST_CASE("Signal Tools") {

            using namespace std::complex_literals;

            const int N = 2;
            const double Wn = 0.5;
            const std::vector<int> x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

            SECTION("Validate Pad") {

                const int ntaps = 3;
                const auto [ext, edge] = detail::_validate_pad<int, PADTYPE::ODD>(x, -1, ntaps);
                REQUIRE(almost_equal<1>(ext, {-9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4,
                    5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}));
                REQUIRE(almost_equal<1>(edge, 9));
                
            }
            SECTION("Companion") {
                const std::vector<double> a = {1., 2., 3., 4., 5.};
                const auto result = detail::companion(a);
                REQUIRE(almost_equal<1>(result, {{-2., -3., -4., -5.},
                                                 { 1.,  0.,  0.,  0.},
                                                 { 0.,  1.,  0.,  0.},
                                                 { 0.,  0.,  1.,  0.}}
                                                ));
                
            }
            SECTION("LFilter ZI") {

                auto ba = butter<BTYPE::LOWPASS, FOUTPUT::BA, double>(2, 0.5);
                const auto result = lfilter_zi(ba.b, ba.a);
                REQUIRE(almost_equal<1000000000>(result, {0.70710678, 0.12132034}));
            }
            SECTION("Raw Filter") {

                auto ba = butter<BTYPE::LOWPASS, FOUTPUT::BA, double>(2, 0.5);
                const auto zi = lfilter_zi(ba.b, ba.a);

                auto [y, zf] = detail::_raw_filter(ba.b, ba.a, {0.0, 1.0, 2.0, 3.0, 4.0}, zi);

                REQUIRE(almost_equal<1000000000>(y, {0.70710678, 0.41421356, 1.05025253, 2.27207794, 3.33452378}));
                REQUIRE(almost_equal<1000000000>(zf, {2.83199846, 0.59945904}));
            }
        }
    }
}
