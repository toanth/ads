
#include "benchmark/benchmark.h"
#include "benchmarks/benchmarks_common.hpp"
#include "include/bitvector/recursive_bitvec.hpp"

using namespace ads;

int main() {
    Span<const U64> randomData = getRandomNumbers(Elem(1) << 28);
    auto bv = EfficientBitvec<>(randomData);
    Index maxVal = bv.numOnes();
    for (Index i = 0;; ++i) {
        Index val = bv.selectOne(randomData[i % randomData.size()] % maxVal);
        benchmark::DoNotOptimize(val);
    }
    return 0;
}
