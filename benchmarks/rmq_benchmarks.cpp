#include "../include/linear_space_rmq.hpp"
#include "../include/nlogn_rmq.hpp"
#include "../include/succinct_rmq.hpp"
#include "benchmarks_common.hpp"

#include <unordered_map>

using namespace ads;

template<typename Rmq>
static void rmqRandomValuesImpl(bm::State& state) {
    Index size = state.max_iterations;
    Rmq rmq(randomArray<Limb>(state.range()));
    std::vector<Index> rmQueries = randomArray<Index>(size * 2, state.range());

    Index i = 0;
    for (auto _ : state) {
        Index l = rmQueries[i++];
        Index u = rmQueries[i++];
        std::tie(l, u) = std::minmax(l, u);
        if (l == u) [[unlikely]] {
            if (l == 0) [[unlikely]] {
                u = 1;
            } else {
                --l;
            }
        }
        auto res = rmq(l, u);
        bm::DoNotOptimize(res);
    }
    setNumBits(state, rmq.allocatedSizeInBits());
    state.SetComplexityN(state.range());
    setGroup(state, 21);
}

static void BM_RmqRandomValuesSucinct(bm::State& state) {
    rmqRandomValuesImpl<SuccinctRMQ<>>(state);
}
static void BM_RmqRandomValuesNLogNSpace(bm::State& state) {
    rmqRandomValuesImpl<NLogNRmq<>>(state);
}
static void BM_RmqRandomValuesLinearSpace(bm::State& state) {
    rmqRandomValuesImpl<LinearSpaceRMQ<>>(state);
}
//
// BENCHMARK(BM_RmqRandomValuesSucinct)->RangeMultiplier(5)->Range(5, maxNumValues / 625)->Complexity(bm::o1);
// BENCHMARK(BM_RmqRandomValuesNLogNSpace)->RangeMultiplier(5)->Range(5, maxNumValues / 625)->Complexity(bm::o1);
// BENCHMARK(BM_RmqRandomValuesLinearSpace)->RangeMultiplier(5)->Range(5, maxNumValues / 625)->Complexity(bm::o1);
