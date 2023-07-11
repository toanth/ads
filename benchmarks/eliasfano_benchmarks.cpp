
#include "../include/elias_fano.hpp"
#include "benchmarks_common.hpp"

#include <unordered_map>

using namespace ads;

// constexpr Index maxNumValues = Index(1) << 27;

static std::vector<Limb> createSortedRandomVector(Index size = maxNumValues) noexcept {
    std::vector<Limb> vec = randomArray(size);
    std::sort(vec.begin(), vec.end());
    return vec;
}

static std::unordered_map<Index, std::vector<Limb>> efVecs /*= createSortedVectors()*/;

static const std::vector<Limb>& getEFNumbers(Index size) noexcept {
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
    std::vector<Limb> numbers = getEFNumbers(state.range());
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
    std::vector<Limb> vec = getEFNumbers(state.range());
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
    std::vector<Limb> predecessorQueries = randomArray<Limb>(size, std::numeric_limits<Limb>::max(), ef.getSmallest());

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
    std::vector<Limb> arr(state.range());
    std::iota(arr.begin(), arr.end(), 42);
    std::vector<Limb> predecessorQueries = randomArray<Limb>(size, state.range() + 42, 42);
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


static void BM_EliasFanoPredecessorAcendingThenLargeRandomQueries(bm::State& state) {
    Index size = state.max_iterations;
    std::vector<Limb> arr(state.range() - 1);
    std::iota(arr.begin(), arr.end(), 0);
    arr.push_back(Limb(-1));
    std::vector<Limb> predecessorQueries = randomArray<Limb>(size, state.range() + state.range() / 10);
    for (Limb& val : predecessorQueries) {
        if (val > state.range()) {
            val -= state.range();
            val = Limb(-1) * (val / (state.range() / 10));
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


static void BM_EliasFanoPredecessorAcendingThenLargeInArray(bm::State& state) {
    Index clusterSize = Index(state.range());
    std::vector<Limb> arr(clusterSize + 1);
    std::iota(arr.begin(), arr.begin() + clusterSize, 0);
    arr.push_back(Limb(-1));
    EliasFano ef(arr);
    std::mt19937 engine(std::random_device{}());
    std::shuffle(arr.begin(), arr.end(), engine);

    Index i = 0;
    for (auto _ : state) {
        auto res = ef.predecessor(arr[i++ % arr.size()]);
        bm::DoNotOptimize(res);
    }
    setNumBits(state, ef.numAllocatedBits());
    state.SetComplexityN(state.range());
    setGroup(state, 11);
}


BENCHMARK(BM_EliasFanoCreation)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
BENCHMARK(BM_EliasFanoPredecessorInArray)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
BENCHMARK(BM_EliasFanoPredecessorRandom)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
BENCHMARK(BM_EliasFanoPredecessorAcending)->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
BENCHMARK(BM_EliasFanoPredecessorAcendingThenLargeRandomQueries())->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
BENCHMARK(BM_EliasFanoPredecessorAcendingThenLargeInArray())->RangeMultiplier(5)->Range(5, maxNumValues)->Complexity(bm::oN);
