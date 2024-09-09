#ifndef SCICPP_SIGNAL_ARRAYTOOLS
#define SCICPP_SIGNAL_ARRAYTOOLS

#include "scicpp/polynomials/polynomial.hpp"
#include "scicpp/core/macros.hpp"
#include "scicpp/core/maths.hpp"
#include "scicpp/core/meta.hpp"
#include "scicpp/core/numeric.hpp"
#include "scicpp/core/utils.hpp"
#include "scicpp/core/range.hpp"
#include "scicpp/core/constants.hpp"
#include "scicpp/core/units/units.hpp"
#include "scicpp/core/print.hpp"
#include <vector>
#include <complex>

namespace scicpp::signal {
    using namespace std;
    namespace detail {
        
        template<typename T>
        auto slice(const std::vector<T>& arr, int start, int stop, int step = 1) {

            std::vector<T> result;
            int n = static_cast<int>(arr.size());

            start = start < 0 ? n + start : start;
            stop = stop < 0 ? n + stop : stop;

            int left = min(start, stop);
            int right = max(start, stop);

            if (left < 0) left = -1;
            if (right > n) right = n;

            while (right > left) {
                if ((step < 0) && (stop < start)) {
                    result.push_back(arr[right]);
                    right += step;
                } else if ((step > 0) && (stop > start)) {
                    result.push_back(arr[left]);
                    left += step;
                } else {
                    break;
                }
            }

            return result;
        }

        template<typename T, typename... Vectors>
        auto concatenate(Vectors&&... vecs) {
            std::vector<T> result;
            (result.insert(
                result.end(),
                std::make_move_iterator(vecs.begin()),
                std::make_move_iterator(vecs.end())
            ), ...);
            return result;
        }

        template <typename T>
        auto eye(size_t N, size_t M = 0, int k = 0) {
            if (M == 0) M = N;
            std::vector<std::vector<T>> result;
            for (size_t i = 0; i < N; ++i) {
                std::vector<T> row = zeros<T>(M);
                if ((i + k >= 0) && (i + k < N)) row[i + k] = 1;
                result.push_back(row);
            }
            return result;
        }
    }

    template<typename T>
    auto axis_slice(const std::vector<T>& arr, int start, int stop, int step = 1) {
        return detail::slice(arr, start, stop, step);
    }

    template<typename T>
    auto odd_ext(const std::vector<T>& x, int n) {
        using namespace scicpp::operators;
        if (n < 1) return x;
        if (n > x.size() - 1) throw std::invalid_argument("Analog filters do not use sampling frequency."); // update
        
        auto left_end = axis_slice<T>(x, 0, 1)[0] * 2;
        auto left_ext = axis_slice<T>(x, n, 0, -1);
        auto right_end = axis_slice<T>(x, -1, x.size())[0] * 2;
        auto right_ext = axis_slice<T>(x, -2, -(n + 2), -1);

        auto ext = detail::concatenate<T>(left_end - left_ext, x, right_end - right_ext);

        return ext;
    }
}

#endif
