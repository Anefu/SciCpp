#ifndef SCICPP_SIGNAL_FILTER
#define SCICPP_SIGNAL_FILTER

#include "scicpp/polynomials/polynomial.hpp"
#include "scicpp/core/macros.hpp"
#include "scicpp/core/maths.hpp"
#include "scicpp/core/meta.hpp"
#include "scicpp/core/numeric.hpp"
#include "scicpp/core/utils.hpp"
#include "scicpp/core/range.hpp"
#include "scicpp/core/constants.hpp"
#include "scicpp/core/units/units.hpp"
#include <vector>
#include <complex>

namespace scicpp::signal {
    enum class BTYPE : int { LOWPASS, HIGHPASS, BANDPASS, BANDSTOP };
    enum class FTYPE : int { BUTTER };
    enum class FOUTPUT : int { BA, ZPK, SOS };

    struct ZPK {
        std::vector<std::complex<double>> z;
        std::vector<std::complex<double>> p;
        double k;
    };

    struct BA {
        std::vector<double> a;
        std::vector<double> b;
    };

    namespace detail {
        auto buttap(int N) {

            using namespace scicpp::units::literals;
            using namespace std::literals::complex_literals;
            // Check N > 0 and is int
            std::vector<std::complex<double>> z;
            std::vector<double> m = arange(-N + 1, N + 1, 2);

            std::vector<std::complex<double>> p = map([=](auto x) {
                return -std::exp(1.0i * x * pi<double> / (2.0 * static_cast<double>(N)));
            }, std::move(m)); // e ^ -(i * pi * Mi) / (2 * N)
            
            double k = 1;

            struct ZPK zpk = { z, p, k };

            return zpk;
        }

        std::size_t _relative_degree(std::size_t z_size, std::size_t p_size) {
            // check diff is < 0
            return p_size - z_size;
        }

        auto bilinear_zpk(const struct ZPK& zpk, double fs) {
            using namespace scicpp::operators;
            std::size_t degree = _relative_degree(zpk.z.size(), zpk.p.size());

            double fs2 = 2.0 * fs;

            std::vector<std::complex<double>> z_z = map([=](auto x) { return ( fs2 + x )/( fs2 - x ); }, zpk.z);
            std::vector<std::complex<double>> p_z = map([=](auto x) { return ( fs2 + x )/( fs2 - x ); }, zpk.p);

            std::vector<double> neg_ones = -1.0 * ones<double>(degree);
            z_z.insert(z_z.end(), neg_ones.begin(), neg_ones.end());

            std::vector<std::complex<double>> fs2_z = fs2 - zpk.z;
            std::vector<std::complex<double>> fs2_p = fs2 - zpk.p;

            double k_z = zpk.k * real(prod(fs2_z)/prod(fs2_p));

            struct ZPK zpk_z = { z_z, p_z, k_z };

            return zpk_z;
        }

        void lp2lp_zpk(struct ZPK& zpk, double wo = 1.0) {
            using namespace scicpp::operators;

            std::size_t degree = _relative_degree(zpk.z.size(), zpk.p.size());

            zpk.z = wo * zpk.z;
            zpk.p = wo * zpk.p;

            zpk.k = zpk.k * std::pow(wo, static_cast<double>(degree));
        }

        auto zpk2tf(const struct ZPK& zpk) {
            using namespace scicpp::polynomial;

            std::vector<double> b;
            std::vector<std::complex<double>> a;

            // check if z is 2D (z is strictly 1D for now)

            // if(zpk.z.size() > 1) {
            //     auto temp = polyfromroots(atleast_1d(zpk.z[0]));
            //     b = zeros<std::complex<double>>(zpk.z.size() + 1);
            //     // auto k = zpk.k * zpk.z.size();
            //     for (int i = 0; i < zpk.z.size(); i++) {
            //         b[i] = zpk.k * polyfromroots(atleast_1d(zpk.z[i])).data()[0]; // use atleast_1d
            //     }
            //     // auto a = map([=](auto x) { polyfromroots(x)[0]; }, zpk.z);
            // } else {
            auto poly_z = polyfromroots(zpk.z);
            b = map([=](auto x) { return zpk.k * std::abs(x); }, poly_z.data());
            a = polyfromroots(zpk.p).data();
            auto real_a = map([=](auto x) { return std::abs(x); }, a);
            // }
            std::reverse(real_a.begin(), real_a.end());
            std::reverse(b.begin(), b.end());

            struct BA ba = {real_a, b};

            return ba;
        }

        template <typename T, FTYPE ftype = FTYPE::BUTTER, BTYPE btype = BTYPE::LOWPASS, FOUTPUT output = FOUTPUT::BA>
        auto iirfilter(int N, double Wn, std::optional<double> rp = std::nullopt,
                       std::optional<double> rs = std::nullopt, bool analog = false, std::optional<double> fs = std::nullopt) {
            double warped;
            double fs_val = 2.0;
            struct ZPK zpk;
            struct BA ba;

            if (fs) {
                if (analog) {
                    throw std::invalid_argument("Analog filters do not use sampling frequency.");
                }
                Wn = 2 * Wn / fs.value();
            }

            if (ftype == FTYPE::BUTTER) {
                zpk = buttap(N);
            } //elif besselap...

            if (!analog) {
                // check 0 < Wn < 1
                fs_val = fs.value_or(2.0);
                warped = 2 * fs_val * scicpp::tan(units::radian<double>(pi<double> * Wn / fs_val));
            } else {
                warped = Wn;
            }

            if (btype == BTYPE::LOWPASS || btype == BTYPE::HIGHPASS) {
                // check size(Wn) == 1
                if (btype == BTYPE::LOWPASS) {
                    lp2lp_zpk(zpk, warped);
                }
            }

            if(!analog) {
               zpk = bilinear_zpk(zpk, fs_val);
            }

            // if (output == ZPK) {
            //     return zpk;
            // } else if (output == BA)
            // {
                ba = zpk2tf(zpk);
            // } else if (output == SOS)
            // {
            //     // return zpk2sos(zpk);
            // }
            return ba;
        }
    }

    template <BTYPE btype = BTYPE::LOWPASS, FOUTPUT output = FOUTPUT::BA, typename T>
    auto butter(int N, double Wn, bool analog = false, std::optional<double> fs = std::nullopt) {
        return detail::iirfilter<T, FTYPE::BUTTER, btype, output>(N, Wn, -1, -1, analog, fs);
    }
} // namespace scicpp::signal

#endif
