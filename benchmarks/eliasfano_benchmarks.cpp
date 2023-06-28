
#include "benchmarks_common.hpp"
#include "include/elias_fano.hpp"

using namespace ads;

constexpr Index maxNumValues = Index(1) << 28;

static std::vector<Elem> createSortedRandomVector() noexcept {
    std::vector<Elem> vec = randomArray(maxNumValues);
    std::sort(vec.begin(), vec.end());
    return vec;
}

static const std::vector<Elem> sorted = createSortedRandomVector();

static Subrange<typename std::vector<Elem>::const_iterator> efNumbers(Index n) noexcept {
    assert(n >= 0 && n <= sorted.size());
    Index offset = (Index(sorted.size()) - n) / 2;
    return {sorted.begin() + offset, sorted.begin() + offset + n};
}

static void BM_EliasFanoCreation(bm::State& state) {
    assert(state.max_iterations <= state.range());
    Subrange r = efNumbers(state.range());
    for (auto _ : state) {
        EliasFano ef(r);
        bm::DoNotOptimize(ef);
        bm::ClobberMemory();
    }
    setNumBits(state, EliasFano(r).numAllocatedBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
}

static void BM_EliasFanoPredecessorInArray(bm::State& state) {
    Subrange r = efNumbers(state.range());
    std::vector<Elem> vec(r.begin(), r.end());
    std::mt19937_64 engine(std::random_device{}());
    std::shuffle(vec.begin(), vec.end(), engine);
    EliasFano ef(r);

    Index i = 0;
    for (auto _ : state) {
        auto res = ef.predecessor(vec[i++ % vec.size()]);
        bm::DoNotOptimize(res);
    }
    setNumBits(state, ef.numAllocatedBits());
    state.SetComplexityN(state.range());
}

static void BM_EliasFanoPredecessorRandom(bm::State& state) {
    Index size = state.max_iterations;
    Subrange r = efNumbers(state.range());
    std::vector<Elem> predecessorQueries = randomArray<Elem>(size, std::numeric_limits<Elem>::max(), *r.begin());
    EliasFano ef(r);

    Index i = 0;
    for (auto _ : state) {
        auto res = ef.predecessor(predecessorQueries[i++]);
        bm::DoNotOptimize(res);
    }
    setNumBits(state, ef.numAllocatedBits());
    state.SetComplexityN(state.range());
}


static void BM_EliasFanoPredecessorAcending(bm::State& state) {
    Index size = state.max_iterations;
    std::vector<Elem> arr(state.range());
    std::iota(arr.begin(), arr.end(), 42);
    std::vector<Elem> predecessorQueries = randomArray<Elem>(size, state.range() + 42, 42);
    EliasFano ef(arr);

    Index i = 0;
    for (auto _ : state) {
        auto res = ef.predecessor(predecessorQueries[i++]);
        bm::DoNotOptimize(res);
    }
    setNumBits(state, ef.numAllocatedBits());
    state.SetComplexityN(state.range());
}



BENCHMARK(BM_EliasFanoCreation)->RangeMultiplier(5)->Range(1, maxNumValues)->Complexity(bm::oN);
BENCHMARK(BM_EliasFanoPredecessorInArray)->RangeMultiplier(5)->Range(1, maxNumValues)->Complexity(bm::oN);
BENCHMARK(BM_EliasFanoPredecessorRandom)->RangeMultiplier(5)->Range(1, maxNumValues)->Complexity(bm::oN);
BENCHMARK(BM_EliasFanoPredecessorAcending)->RangeMultiplier(5)->Range(1, maxNumValues)->Complexity(bm::oN);
