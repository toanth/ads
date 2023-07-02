#ifndef ADS_BENCHMARKS_COMMON_HPP
#define ADS_BENCHMARKS_COMMON_HPP

#include "benchmark/benchmark.h"
#include "include/common.hpp"
#include <unordered_map>

namespace bm = benchmark;

namespace ads {


constexpr Index maxNumValues = Index(1) << 27;

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

static std::unordered_map<Index, std::vector<U64>> vecs;

// Computing fresh random numbers for each benchmark invocation can take a long time, so only do the work once (per
// file). Usually, only benchmarks from one file are run, so the fact everything is `static` shouldn't matter.
static Span<const U64> getRandomNumbers(Index size) noexcept {
    auto it = vecs.find(size);
    if (it == vecs.end()) {
        it = vecs.try_emplace(size, randomArray(size)).first;
    }
    return it->second;
}



static void divideByNInPlot(bm::State& state) {
    state.counters["perN"] += 1;
}

static void subtractNFromBitCount(bm::State& state) {
    state.counters["perN"] += 2;
}

static void setNumBits(bm::State& state, Index numBits) {
    state.counters["bits"] = double(numBits);
}

static void setGroup(bm::State& state, double groupIdx) {
    state.counters["group"] = groupIdx;
}


} // namespace ads

#endif // ADS_BENCHMARKS_COMMON_HPP
