#include "benchmark/benchmark.h"
#include "include/bitvector.hpp"

namespace bm = benchmark;
using namespace ads;


static void BM_BitvecAllocation(bm::State& state) {
    for (auto _: state) {
        Bitvector<> bv(state.range(0));
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
}

static void BM_BitvecFillZero(bm::State& state) {
    for (auto _: state) {
        Bitvector<> bv(state.range(0), 0);
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
}

static void BM_BitvecFillOnes(bm::State& state) {
    for (auto _: state) {
        Bitvector<> bv(state.range(0), Elem(-1));
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
}

constexpr Elem alternating = 0xaaaa'aaaa'aaaa'aaaaull;

static void BM_BitvecAlternatingOneZerosRankLast(bm::State& state) {
    Bitvector<> bv(state.range(0), alternating);
    for (auto _: state) {
        Index i = bv.rankZero(state.range(0) - 1);
        bm::DoNotOptimize(i);
    }
}

static void BM_BitvecAlternatingOneZerosSelectLast(bm::State& state) {
    Bitvector<> bv(state.range(0), alternating);
    for (auto _: state) {
        Index i = bv.selectZero(state.range(0) / 2);
        bm::DoNotOptimize(i);
    }
}


BENCHMARK(BM_BitvecAllocation)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecFillZero)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecFillOnes)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecAlternatingOneZerosRankLast)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecAlternatingOneZerosSelectLast)->Range(1, Elem(1) << 34);
