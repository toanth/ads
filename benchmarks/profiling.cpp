
#include "benchmark/benchmark.h"
#include "benchmarks/benchmarks_common.hpp"
#include "include/bitvector/recursive_bitvec.hpp"
#include "include/elias_fano.hpp"

#include <iostream>

using namespace ads;

int bitvec() {
    Span<const U64> randomData = getRandomNumbers(U64(1) << 28);
    auto bv = EfficientBitvec<>(randomData);
    Index maxVal = bv.numOnes();
    for (Index i = 0;; ++i) {
        Index val = bv.selectOne(randomData[i % randomData.size()] % maxVal);
        benchmark::DoNotOptimize(val);
    }
    return 0;
}

int eliasFano() {
    std::cout << "Profiling elias fano.....set up" << std::endl;
    auto unsortedRandomData = getRandomNumbers(U64(1) << 28);
    std::vector<U64> randomData(unsortedRandomData.begin(), unsortedRandomData.end());
    randomData.push_back(0);
    std::sort(randomData.begin(), randomData.end());
    EliasFano ef(randomData);
    std::cout << "Starting infinite loop" << std::endl;
    for (Index i = 0;; ++i) {
        Index val = ef.predecessor(randomData[i % randomData.size()]);
        benchmark::DoNotOptimize(val);
    }
    return 0;
}

int main() {
    eliasFano();
}
