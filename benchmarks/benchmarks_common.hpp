#ifndef ADS_BENCHMARKS_COMMON_HPP
#define ADS_BENCHMARKS_COMMON_HPP

#include "benchmark/benchmark.h"
#include "include/common.hpp"

namespace bm = benchmark;

namespace ads {

template<typename T = Elem>
static std::vector<T> randomArray(Index size, T maxVal = std::numeric_limits<T>::max(), T minVal = 0) {
    std::vector<T> vec(size);
    std::uniform_int_distribution<T> dist(minVal, maxVal);
    std::random_device rd;
    std::mt19937_64 engine(rd());
    for (auto& v : vec) {
        v = dist(engine);
    }
    return vec;
}


static void divideByNInPlot(bm::State& state) {
    state.counters["perN"] = 1;
}

static void subtractNAndDivideByNInPlot(bm::State& state) {
    state.counters["perN"] = 2;
}

static void setNumBits(bm::State& state, Index numBits) {
    state.counters["bits"] = double(numBits);
}

static void setGroup(bm::State& state, double groupIdx) {
    state.counters["group"] = groupIdx;
}


} // namespace ads

#endif // ADS_BENCHMARKS_COMMON_HPP
