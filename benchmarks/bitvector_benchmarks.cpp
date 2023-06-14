#include "benchmark/benchmark.h"
#include "include/bitvector.hpp"

namespace bm = benchmark;
using namespace ads;


static void BM_BitvecAllocation(bm::State& state) {
    for (auto _ : state) {
        Bitvector<> bv(state.range(0));
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    state.counters["bits"] = double(Bitvector<>(state.range(0)).numAllocatedBits());
}

static void BM_BitvecFillZero(bm::State& state) {
    for (auto _ : state) {
        Bitvector<> bv(state.range(0), 0);
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    state.counters["bits"] = double(Bitvector<>(state.range(0)).numAllocatedBits());
}

static void BM_BitvecFillOnes(bm::State& state) {
    for (auto _ : state) {
        Bitvector<> bv(state.range(0), Elem(-1));
        bm::DoNotOptimize(bv);
        bm::ClobberMemory();
    }
    state.counters["bits"] = double(Bitvector<>(state.range(0)).numAllocatedBits());
}

constexpr Elem alternating = 0xaaaa'aaaa'aaaa'aaaaull;

static void BM_BitvecAlternatingOnesZerosRankLast(bm::State& state) {
    Bitvector<> bv(state.range(0), alternating);
    for (auto _ : state) {
        Index i = bv.rankZero(state.range(0) - 1);
        bm::DoNotOptimize(i);
    }
    state.counters["bits"] = double(bv.numAllocatedBits());
}

static void BM_BitvecAlternatingOnesZerosSelectLast(bm::State& state) {
    Bitvector<> bv(state.range(0), alternating);
    for (auto _ : state) {
        Index i = bv.selectZero(state.range(0) / 2);
        bm::DoNotOptimize(i);
    }
    state.counters["bits"] = double(bv.numAllocatedBits());
}

static void BM_BitvecAlternatingOnesZerosSelectOneThird(bm::State& state) {
    Bitvector<> bv(state.range(0), alternating);
    for (auto _ : state) {
        Index i = bv.selectOne(state.range(0) / 3);
        bm::DoNotOptimize(i);
    }
    state.counters["bits"] = double(bv.numAllocatedBits());
}

static void BM_BitvecOnesThenZerosSelectFirstZero(bm::State& state) {
    Bitvector<> bv(state.range(0));
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
    state.counters["bits"] = double(bv.numAllocatedBits());
}

static void BM_Bitvec4OnesThenZerosThenOnesSelectFifthOne(bm::State& state) {
    Bitvector<> bv(state.range(0) + 64);
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
    state.counters["bits"] = double(bv.numAllocatedBits());
}

BENCHMARK(BM_BitvecAllocation)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecFillZero)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecFillOnes)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecAlternatingOnesZerosRankLast)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecAlternatingOnesZerosSelectLast)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecAlternatingOnesZerosSelectOneThird)->Range(1, Elem(1) << 34);
BENCHMARK(BM_BitvecOnesThenZerosSelectFirstZero)->Range(1, Elem(1) << 34);
BENCHMARK(BM_Bitvec4OnesThenZerosThenOnesSelectFifthOne)->Range(1, Elem(1) << 34);
