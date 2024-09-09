#include "filter_design.hpp"

#include "scicpp/core/equal.hpp"
#include "scicpp/core/numeric.hpp"
#include "scicpp/core/print.hpp"

namespace scicpp {
    namespace signal {
        TEST_CASE("Butterworth Filter") {

            using namespace std::complex_literals;

            const int N = 2;
            const double Wn = 0.5;

            SECTION("Buttap") {

                const auto result = detail::buttap(N);

                REQUIRE(result.z.empty());
                REQUIRE(almost_equal<100000000>(result.p, {-0.70710678 + 0.70710678i, -0.70710678 - 0.70710678i}));
                REQUIRE(almost_equal<1>(result.k, 1.0));
            }

            SECTION("LP2LP ZPK") {

                auto zpk = detail::buttap(N);
                double warped = 2 * 2.0 * scicpp::tan(units::radian<double>(pi<double> * Wn / 2.0));

                // const auto result = detail::lp2lp_zpk(zpk, warped);
                detail::lp2lp_zpk(zpk, warped);

                REQUIRE(zpk.z.empty());
                REQUIRE(almost_equal<100000000>(zpk.p, {-2.82842712+2.82842712i, -2.82842712-2.82842712i}));
                REQUIRE(almost_equal<1>(zpk.k, 15.999999999999996));
            }

            SECTION("Bilinear ZPK") {

                auto zpk_b = detail::buttap(N);
                double warped = 2 * 2.0 * scicpp::tan(units::radian<double>(pi<double> * Wn / 2.0));

                // const auto zpk_l = detail::lp2lp_zpk(zpk_b, warped);
                detail::lp2lp_zpk(zpk_b, warped);

                const auto result = detail::bilinear_zpk(zpk_b, 2.0);

                scicpp::print(zpk_b.p);

                REQUIRE(almost_equal<1>(result.z, {-1., -1.}));
                REQUIRE(almost_equal<100000000>(result.p, {9.20184753e-17+0.41421356i, 9.20184753e-17-0.41421356i}));
                REQUIRE(almost_equal<100000000>(result.k, 0.2928932188134525));
            }

            SECTION("Lowpass Filter with N = 2, Wn = 0.5") {
                const bool analog = false;
                const auto result = butter<BTYPE::LOWPASS, FOUTPUT::BA, double>(N, Wn);

                REQUIRE(almost_equal<100000000>(result.b, {0.29289322, 0.58578644, 0.29289322}));
                REQUIRE(almost_equal<100000000>(result.a, {1.00000000e+00, -1.84036951e-16,  1.71572875e-01}));
            }

            // SECTION("Highpass Filter with Analog True") {
            //     const int N = 2;
            //     const double Wn = 1.0;
            //     const bool analog = true;
            //     const double fs = -1;
            //     const auto result = butter<HIGHPASS, BA, double>(N, Wn, analog, fs);
                
            //     REQUIRE(almost_equal(result.p, {1.0i, 2.0i})); // Replace with actual expected values
            //     REQUIRE(result.k == Approx(1.0)); // Replace with actual expected value
            // }

            // SECTION("Bandpass Filter with Sampling Frequency") {
            //     const int N = 4;
            //     const double Wn = 0.25;
            //     const bool analog = false;
            //     const double fs = 1000;
            //     const auto result = butter<BANDPASS, BA, double>(N, Wn, analog, fs);
                
            //     REQUIRE(almost_equal(result.p, {1.0i, 2.0i, 3.0i, 4.0i})); // Replace with actual expected values
            //     REQUIRE(result.k == Approx(1.0)); // Replace with actual expected value
            // }

            // SECTION("Bandstop Filter with Analog False and Sampling Frequency") {
            //     const int N = 3;
            //     const double Wn = 0.3;
            //     const bool analog = false;
            //     const double fs = 2000;
            //     const auto result = butter<BANDSTOP, BA, double>(N, Wn, analog, fs);
                
            //     REQUIRE(almost_equal(result.p, {1.0i, 2.0i, 3.0i})); // Replace with actual expected values
            //     REQUIRE(result.k == Approx(1.0)); // Replace with actual expected value
            // }
        }
    }
}