
#include "benchmarks_common.hpp"
#include "include/elias_fano.hpp"

using namespace ads;

constexpr Index maxNumValues = Index(1) << 31;


static void BM_EliasFanoCreation(bm::State& state) {
    std::vector<Elem> vec = randomArray(state.range());
    std::sort(vec.begin(), vec.end());
    for (auto _ : state) {
        EliasFano ef(vec);
        bm::DoNotOptimize(ef);
        bm::ClobberMemory();
    }
    setNumBits(state, EliasFano(vec).numAllocatedBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
}

BENCHMARK(BM_EliasFanoCreation)->RangeMultiplier(5)->Range(1, maxNumValues)->Complexity(bm::oN);
