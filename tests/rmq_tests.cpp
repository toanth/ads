
//#include "../include/cartesian_tree_rmq.hpp"
#include "../include/naive_rmq.hpp"
#include "../include/nlogn_rmq.hpp"
#include "gtest/gtest.h"

using namespace ads;

template<typename RMQ>
void testSmallRmq() {
    std::vector<Elem> values{12, 0, 4, 2, 100, 101};
    RMQ rmq(values);
    SimpleRMQ<> reference(values);
    for (Index i = 0; i < values.size(); ++i) {
        for (Index j = i + 1; j <= values.size(); ++j) {
            ASSERT_EQ(rmq(i, j), reference(i, j)) << i << " " << j;
        }
    }
}


template<typename RMQ>
void testRandomRmq() {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist;
    std::vector<Elem> values(dist(engine) % (1 << 10));
    for (auto& val: values) {
        val = dist(engine);
    }
    RMQ rmq(values);
    SimpleRMQ<> reference(values);
    std::uniform_int_distribution<Index> idxDist(0, Index(values.size()));
    for (Index i = 0; i < 100'000; ++i) {
        std::pair<Index, Index> range = std::minmax(idxDist(engine), idxDist(engine));
        if (range.first == range.second || range.first == values.size()) continue;
        ASSERT_EQ(rmq(range.first, range.second), reference(range.first, range.second)) << range.first << " " << range.second;
    }
}
TEST(NaiveRmq, Small) {
    testSmallRmq<NaiveRMQ<>>();
}

TEST(NaiveRmq, Random) {
    testRandomRmq<NaiveRMQ<>>();
}

TEST(NlognRMQ, Small) {
    testSmallRmq<NlognRMQ<>>();
}

TEST(NlognRMQ, Random) {
    testRandomRmq<NlognRMQ<>>();
}

//TEST(CartesianTreeRMQ, Small) {
//    testSmallRmq<CartesianTreeRMQ<>>();
//}
