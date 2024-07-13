#ifndef SCICPP_CORE_APPEND
#define SCICPP_CORE_APPEND

#include <vector>

namespace scicpp {
    template <typename T>
    auto append(std::vector<T>& first, std::vector<T>& second) {
        return first.insert(first_vector.end(), second_vector.begin(), second.end());
    }
}

#endif

