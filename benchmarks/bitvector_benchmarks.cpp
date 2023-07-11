#include "../include/bitvector/recursive_bitvec.hpp"
#include "benchmarks_common.hpp"


using namespace ads;

// using RankBitvector = EfficientRankBitvec<>;
using RankBitvector = EfficientBitvec<>; // also use almost the same bitvector for fairer comparisons
using SelectBitvector = EfficientSelectBitvec<>;
using DefaultBitvector = SelectBitvector;


constexpr Index maxNumBits = maxNumValues * 64 * 2;

// Only used internally, so a macro is fine here.
// Also, the fact that random numbers aren't perfectly uniformly distributed doesn't matter
#define ADS_GET_RANDVAL(maxVal) Index(randomQueries[i++] % maxVal)


[[nodiscard]] Span<const U64> getRandomQueries(const bm::State& state) {
    return getRandomNumbers(state.max_iterations);
}

static void BM_BitvecAllocation(bm::State& state) {
    for (auto _ : state) {
        DefaultBitvector bv = DefaultBitvector::uninitializedForSize(state.range());
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    setNumBits(state, DefaultBitvector::uninitializedForSize(state.range()).allocatedSizeInBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
    subtractNFromBitCount(state);
    setGroup(state, 1);
}

static void BM_BitvecFillZero(bm::State& state) {
    for (auto _ : state) {
        DefaultBitvector bv(state.range(), Limb(0));
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    setNumBits(state, DefaultBitvector(state.range(), 0).allocatedSizeInBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
    subtractNFromBitCount(state);
    setGroup(state, 1);
}

static void BM_BitvecFillOnes(bm::State& state) {
    for (auto _ : state) {
        DefaultBitvector bv(state.range(), Limb(-1));
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    setNumBits(state, DefaultBitvector(state.range(), Limb(-1)).allocatedSizeInBits());
    state.SetComplexityN(state.range());
    divideByNInPlot(state);
    subtractNFromBitCount(state);
    setGroup(state, 1);
}

constexpr Limb alternating = 0xaaaa'aaaa'aaaa'aaaaull;

static void BM_BitvecAlternatingOnesZerosRankTwoThird(bm::State& state) {
    RankBitvector bv(state.range(), alternating);
    for (auto _ : state) {
        Index i = bv.rankZeroUnchecked(state.range() * 2 / 3);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.1);
}

static void BM_BitvecAlternatingOnesZerosRankRandom(bm::State& state) {
    Span<const U64> randomQueries = getRandomQueries(state);
    RankBitvector bv(state.range(), alternating);
    Index i = 0;
    for (auto _ : state) {
        Index val = bv.rankZeroUnchecked(ADS_GET_RANDVAL(state.range()));
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.1);
}

static void BM_BitvecRandomRank(bm::State& state) {
    Span<const U64> randomQueries = getRandomQueries(state);
    assert(state.range() <= maxNumBits);
    RankBitvector bv(getRandomNumbers(roundUpDiv(state.range(), 64)), state.range());
    Index i = 0;
    for (auto _ : state) {
        Index val = bv.rankOneUnchecked(ADS_GET_RANDVAL(state.range()));
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.1);
}


// ** Select **

static void BM_BitvecAlternatingOnesZerosSelectOneThird(bm::State& state) {
    SelectBitvector bv(state.range(), alternating);
    for (auto _ : state) {
        Index i = bv.selectOne(state.range() / 3);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.2);
}

static void BM_BitvecAlternatingOnesZerosSelectRandom(bm::State& state) {
    Span<const U64> randomQueries = getRandomQueries(state);
    SelectBitvector bv(state.range(), alternating);
    Index maxRank = bv.numZeros();
    assert(maxRank > 0);
    Index i = 0;
    for (auto _ : state) {
        Index val = bv.selectZero(ADS_GET_RANDVAL(maxRank));
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.2);
}

static void BM_BitvecOnesThenZerosSelectFirstZero(bm::State& state) {
    SelectBitvector bv = SelectBitvector::uninitializedForSize(state.range());
    for (Index i = 0; i < bv.numLimbs() / 2; ++i) {
        bv.setLimb(i, Limb(-1));
    }
    for (Index i = bv.numLimbs() / 2; i < bv.numLimbs(); ++i) {
        bv.setLimb(i, 0);
    }
    bv.buildMetadata();
    for (auto _ : state) {
        Index i = bv.selectZero(0);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.2);
}

static void BM_BitvecOnesThenZerosSelectRandom(bm::State& state) {
    Span<const U64> randomQueries = getRandomQueries(state);
    SelectBitvector bv = SelectBitvector::uninitializedForSize(state.range());
    for (Index i = 0; i < bv.numLimbs() / 2; ++i) {
        bv.setLimb(i, Limb(-1));
    }
    for (Index i = bv.numLimbs() / 2; i < bv.numLimbs(); ++i) {
        bv.setLimb(i, 0);
    }
    bv.buildMetadata();
    Index maxRank = bv.numZeros();
    Index i = 0;
    for (auto _ : state) {
        Index val = bv.selectZero(ADS_GET_RANDVAL(maxRank));
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.2);
}

static void BM_BitvecFirstLastBitOneOthersZeroSelectSecondOne(bm::State& state) {
    SelectBitvector bv(state.range(), 0);
    ADS_ASSUME(state.range() >= 2);
    bv.setBit(0);
    bv.setBit(bv.size() - 1);
    bv.buildMetadata();
    if (!bv.getBit(0)) throw std::logic_error("internal error");
    for (auto _ : state) {
        Index i = bv.selectOne(1);
        bm::DoNotOptimize(i);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.2);
}


static void BM_BitvecFirstLastBitOneOthersZeroSelectRandom(bm::State& state) {
    Span<const U64> randomQueries = getRandomQueries(state);
    SelectBitvector bv(state.range(), 0);
    ADS_ASSUME(state.range() >= 2);
    bv.setBit(0);
    bv.setBit(bv.size() - 1);
    bv.buildMetadata();
    if (!bv.getBit(0)) throw std::logic_error("internal error");
    Index i = 0;
    for (auto _ : state) {
        Index val = bv.selectOne(ADS_GET_RANDVAL(2));
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.2);
}

static void BM_BitvecSqrtNRandomOnesSelectRandom(bm::State& state) {
    Span<const U64> randomQueries = getRandomQueries(state);
    SelectBitvector bv(state.range(), 0);
    Index numOnes = std::clamp(Index(1), Index(std::sqrt(state.range())), randomQueries.size());
    ADS_ASSUME(bv.size() > 0);
    for (Index i = 0; i < numOnes; ++i) {
        bv.setBit(randomQueries[i] % bv.size()); // it's ok to set the same bit twice, shouldn't happen too often
    }
    ADS_ASSUME(numOnes >= bv.numOnes());
    numOnes = bv.numOnes();
    bv.buildMetadata();
    Index i = 0;
    for (auto _ : state) {
        Index val = bv.selectOne(ADS_GET_RANDVAL(2));
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.2);
}

static void BM_BitvecSelectRandom(bm::State& state) {
    Span<const U64> randomQueries = getRandomQueries(state);
    SelectBitvector bv(getRandomNumbers(roundUpDiv(state.range(), 64)), state.range());
    if (bv.numOnes() == 0) {
        bv.setBit(0);
        bv.buildMetadata();
    }
    Index i = 0;
    Index maxRank = bv.numOnes();
    for (auto _ : state) {
        Index val = bv.selectOne(ADS_GET_RANDVAL(maxRank));
        bm::DoNotOptimize(val);
    }
    setNumBits(state, bv.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    subtractNFromBitCount(state);
    setGroup(state, 2.2);
}

//
// BENCHMARK(BM_BitvecAllocation)->Range(8, Limb(1) << 34)->Complexity(bm::oN);
// BENCHMARK(BM_BitvecFillZero)->Range(8, Limb(1) << 34)->Complexity(bm::oN);
//// BENCHMARK(BM_BitvecFillOnes)->Range(8, Limb(1) << 34)->Complexity(bm::oN);
//
//
// BENCHMARK(BM_BitvecAlternatingOnesZerosRankTwoThird)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::o1);
// BENCHMARK(BM_BitvecAlternatingOnesZerosRankRandom)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::o1);
// BENCHMARK(BM_BitvecRandomRank)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::o1);


// BENCHMARK(BM_BitvecAlternatingOnesZerosSelectOneThird)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::o1);
// BENCHMARK(BM_BitvecAlternatingOnesZerosSelectRandom)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::o1);

// BENCHMARK(BM_BitvecOnesThenZerosSelectFirstZero)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::o1);
BENCHMARK(BM_BitvecOnesThenZerosSelectRandom)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::o1);

// BENCHMARK(BM_BitvecFirstLastBitOneOthersZeroSelectSecondOne)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::oLogN);
BENCHMARK(BM_BitvecFirstLastBitOneOthersZeroSelectRandom)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::oLogN);

BENCHMARK(BM_BitvecSqrtNRandomOnesSelectRandom)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::oLogN);

BENCHMARK(BM_BitvecSelectRandom)->RangeMultiplier(5)->Range(5, maxNumBits)->Complexity(bm::o1);
