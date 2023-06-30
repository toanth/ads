#include "benchmarks_common.hpp"
#include "include/bitvector.hpp"

using namespace ads;

constexpr Index maxNumBits = Index(1) << 34;

static std::vector<Elem> randomBitvecValues = randomArray(roundUpDiv(maxNumBits, 64));

static void BM_BitvecAllocation(bm::State& state) {
    for (auto _ : state) {
        Bitvector<> bv(state.range());
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    setNumBits(state, Bitvector<>(state.range()).numAllocatedBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
}

static void BM_BitvecFillZero(bm::State& state) {
    for (auto _ : state) {
        Bitvector<> bv(state.range(), Elem(0));
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    setNumBits(state, Bitvector<>(state.range()).numAllocatedBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
}

static void BM_BitvecFillOnes(bm::State& state) {
    for (auto _ : state) {
        Bitvector<> bv(state.range(), Elem(-1));
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    setNumBits(state, Bitvector<>(state.range()).numAllocatedBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
}

constexpr Elem alternating = 0xaaaa'aaaa'aaaa'aaaaull;

static void BM_BitvecAlternatingOnesZerosRankLast(bm::State& state) {
    Bitvector<> bv(state.range(), alternating);
    for (auto _ : state) {
        Index i = bv.rankZeroUnchecked(state.range() - 1);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.numAllocatedBits());
    state.SetComplexityN(state.range());
}

static void BM_BitvecAlternatingOnesZerosRankTwoThird(bm::State& state) {
    Bitvector<> bv(state.range(), alternating);
    for (auto _ : state) {
        Index i = bv.rankZeroUnchecked(state.range() * 2 / 3);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.numAllocatedBits());
    state.SetComplexityN(state.range());
}

static void BM_BitvecRandomRank(bm::State& state) {
    std::vector<Index> randomValues = randomArray<Index>(state.max_iterations, state.range() - 1);
    assert(state.range() <= maxNumBits);
    Bitvector<> bv(randomBitvecValues, state.range());
    Index i = 0;
    for (auto _ : state) {
        Index val = bv.rankOneUnchecked(randomValues[i++]);
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.numAllocatedBits());
    state.SetComplexityN(state.range());
}

static void BM_BitvecAlternatingOnesZerosSelectLast(bm::State& state) {
    Bitvector<> bv(state.range(), alternating);
    for (auto _ : state) {
        Index i = bv.selectZero(state.range() / 2);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.numAllocatedBits());
    state.SetComplexityN(state.range());
}

static void BM_BitvecAlternatingOnesZerosSelectOneThird(bm::State& state) {
    Bitvector<> bv(state.range(), alternating);
    for (auto _ : state) {
        Index i = bv.selectOne(state.range() / 3);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.numAllocatedBits());
    state.SetComplexityN(state.range());
}

static void BM_BitvecOnesThenZerosSelectFirstZero(bm::State& state) {
    Bitvector<> bv(state.range());
    for (Index i = 0; i < bv.sizeInElems() / 2; ++i) {
        bv.setElem(i, Elem(-1));
    }
    for (Index i = bv.sizeInElems() / 2; i < bv.sizeInElems(); ++i) {
        bv.setElem(i, 0);
    }
    bv.buildMetadata();
    for (auto _ : state) {
        Index i = bv.selectZero(0);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.numAllocatedBits());
    state.SetComplexityN(state.range());
}

static void BM_Bitvec4OnesThenZerosThenOnesSelectFifthOne(bm::State& state) {
    Bitvector<> bv(state.range() + 64);
    bv.setElem(0, 0xf);
    if (!bv.getBit(0)) throw std::logic_error("internal error");
    for (Index i = 1; i < bv.sizeInElems() / 2; ++i) {
        bv.setElem(i, 0);
    }
    for (Index i = bv.sizeInElems() / 2; i < bv.sizeInElems(); ++i) {
        bv.setElem(i, Elem(-1));
    }
    bv.buildMetadata();
    for (auto _ : state) {
        Index i = bv.selectOne(5);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.numAllocatedBits());
    state.SetComplexityN(state.range());
}

static void BM_BitvecSelectRandom(bm::State& state) {
    Bitvector<> bv(randomBitvecValues, state.range());
    Index maxRank;
    while (true) {
        maxRank = bv.rankOne(bv.sizeInBits() - 1);
        if (!bv.getBit(bv.sizeInBits() - 1)) {
            --maxRank;
        }
        if (maxRank >= 0) {
            break;
        }
        bv.setElem(0, bv.getElem(0) + 1);
    }
    Index size = state.max_iterations;
    std::vector<Index> selectQueries = randomArray<Index>(size, maxRank);
    Index i = 0;
    for (auto _ : state) {
        Index val = bv.selectOne(selectQueries[i++]);
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.numAllocatedBits());
    state.SetComplexityN(state.range());
}

// BENCHMARK(BM_BitvecAllocation)->Range(1, Elem(1) << 34)->Complexity(bm::oN);
// BENCHMARK(BM_BitvecFillZero)->Range(1, Elem(1) << 34)->Complexity(bm::oN);
// BENCHMARK(BM_BitvecFillOnes)->Range(1, Elem(1) << 34)->Complexity(bm::oN);

// BENCHMARK(BM_BitvecAlternatingOnesZerosRankLast)->Range(1, maxNumBits)->Complexity(bm::o1);
// BENCHMARK(BM_BitvecAlternatingOnesZerosRankTwoThird)->RangeMultiplier(5)->Range(1, maxNumBits)->Complexity(bm::o1);
// BENCHMARK(BM_BitvecRandomRank)->RangeMultiplier(5)->Range(1, maxNumBits)->Complexity(bm::o1);

// BENCHMARK(BM_BitvecAlternatingOnesZerosSelectLast)->RangeMultiplier(5)->Range(1, maxNumBits)->Complexity(bm::oLogN);
// BENCHMARK(BM_BitvecAlternatingOnesZerosSelectOneThird)->RangeMultiplier(5)->Range(1, maxNumBits)->Complexity(bm::oLogN);
// BENCHMARK(BM_BitvecOnesThenZerosSelectFirstZero)->RangeMultiplier(5)->Range(1, maxNumBits)->Complexity(bm::oLogN);
// BENCHMARK(BM_Bitvec4OnesThenZerosThenOnesSelectFifthOne)->RangeMultiplier(5)->Range(1, maxNumBits)->Complexity(bm::oLogN);
// BENCHMARK(BM_BitvecSelectRandom)->RangeMultiplier(5)->Range(1, maxNumBits)->Complexity(bm::oLogN);
