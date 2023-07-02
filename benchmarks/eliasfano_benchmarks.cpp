
#include "benchmarks_common.hpp"
#include "include/elias_fano.hpp"

#include <unordered_map>

using namespace ads;

// constexpr Index maxNumValues = Index(1) << 27;

static std::vector<Elem> createSortedRandomVector(Index size = maxNumValues) noexcept {
    std::vector<Elem> vec = randomArray(size);
    std::sort(vec.begin(), vec.end());
    return vec;
}

static std::unordered_map<Index, std::vector<Elem>> efVecs /*= createSortedVectors()*/;

static const std::vector<Elem>& getEFNumbers(Index size) noexcept {
    auto it = efVecs.find(size);
    if (it == efVecs.end()) {
        it = efVecs.try_emplace(size, createSortedRandomVector(size)).first;
    }
    return it->second;
}

static EliasFano<> getEF(Index size) noexcept {
    return EliasFano(getEFNumbers(size));
}



static void BM_EliasFanoCreation(bm::State& state) {
    std::vector<Elem> numbers = getEFNumbers(state.range());
    for (auto _ : state) {
        EliasFano ef(numbers);
        bm::DoNotOptimize(ef);
        bm::ClobberMemory();
    }
    setNumBits(state, EliasFano(numbers).numAllocatedBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
    setGroup(state, 10);
}

static void BM_EliasFanoPredecessorInArray(bm::State& state) {
    std::vector<Elem> vec = getEFNumbers(state.range());
    EliasFano<> ef(vec);
    std::mt19937_64 engine(std::random_device{}());
    std::shuffle(vec.begin(), vec.end(), engine);

    Index i = 0;
    for (auto _ : state) {
        auto res = ef.predecessor(vec[i++ % vec.size()]);
        bm::DoNotOptimize(res);
    }
    setNumBits(state, ef.numAllocatedBits());
    state.SetComplexityN(state.range());
    setGroup(state, 11);
}

static void BM_EliasFanoPredecessorRandom(bm::State& state) {
    Index size = state.max_iterations;
    EliasFano<> ef = getEF(state.range());
    std::vector<Elem> predecessorQueries = randomArray<Elem>(size, std::numeric_limits<Elem>::max(), ef.getSmallest());

    Index i = 0;
    for (auto _ : state) {
        auto res = ef.predecessor(predecessorQueries[i++]);
        bm::DoNotOptimize(res);
    }
    setNumBits(state, ef.numAllocatedBits());
    state.SetComplexityN(state.range());
    setGroup(state, 11);
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
    setGroup(state, 11);
}


static void BM_EliasFanoPredecessorAcendingThenLarge(bm::State& state) {
    Index size = state.max_iterations;
    std::vector<Elem> arr(state.range() - 1);
    std::iota(arr.begin(), arr.end(), 0);
    arr.push_back(Elem(-1));
    std::vector<Elem> predecessorQueries = randomArray<Elem>(size, state.range() + state.range() / 10);
    for (Elem& val : predecessorQueries) {
        if (val > state.range()) {
            val -= state.range();
            val = Elem(-1) * (val / (state.range() / 10));
        }
    }
    EliasFano ef(arr);

    Index i = 0;
    for (auto _ : state) {
        auto res = ef.predecessor(predecessorQueries[i++]);
        bm::DoNotOptimize(res);
    }
    setNumBits(state, ef.numAllocatedBits());
    state.SetComplexityN(state.range());
    setGroup(state, 11);
}


// BENCHMARK(BM_EliasFanoCreation)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
// BENCHMARK(BM_EliasFanoPredecessorInArray)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
// BENCHMARK(BM_EliasFanoPredecessorRandom)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
// BENCHMARK(BM_EliasFanoPredecessorAcending)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
// BENCHMARK(BM_EliasFanoPredecessorAcendingThenLarge)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
