
//#include "../include/cartesian_tree_rmq.hpp"
#include "../include/linear_space_rmq.hpp"
#include "../include/naive_rmq.hpp"
#include "../include/nlogn_rmq.hpp"
#include "gtest/gtest.h"

using namespace ads;


template<ADS_RMQ_CONCEPT RmqType>
class AllRmqsTest : public ::testing::Test {
};

template<ADS_RMQ_CONCEPT RmqType>
class UsefulRmqsTest : public ::testing::Test {
};

using AllRmqTypes = ::testing::Types<NaiveRMQ<>, NlognRMQ<>, LinearSpaceRMQ<>>;
TYPED_TEST_SUITE(AllRmqsTest, AllRmqTypes);

using UsefulRmqTypes = ::testing::Types<NlognRMQ<>, LinearSpaceRMQ<>>;
TYPED_TEST_SUITE(UsefulRmqsTest, UsefulRmqTypes);


void fillWithRandomValues(std::vector<Elem>& v) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist;
    for (Elem& value: v) {
        value = dist(engine);
    }
}

std::vector<Elem> randomValues(Index size) {
    std::vector<Elem> res(size);
    fillWithRandomValues(res);
    return res;
}

template<ADS_RMQ_CONCEPT TestRmq, ADS_RMQ_CONCEPT ReferenceRmq>
void testQuery(const TestRmq& testRmq, const ReferenceRmq& referenceRmq, Index lower, Index upper) {
    Index testeeAnswer = testRmq(lower, upper);
    Index referenceAnswer = referenceRmq(lower, upper);
    ASSERT_EQ(referenceRmq[testeeAnswer], referenceRmq[referenceAnswer]) << TestRmq::name
                                                                         << ", answers " << testeeAnswer << " "
                                                                         << referenceAnswer << ", interval " << lower
                                                                         << " " << upper;
}

template<typename RMQ>
void testSmallRmq(std::vector<Elem> values) {
    RMQ rmq(values);
    SimpleRMQ<> reference(values);
    for (Index i = 0; i < values.size(); ++i) {
        for (Index j = i + 1; j <= values.size(); ++j) {
            testQuery(rmq, reference, i, j);
        }
    }
}

template<typename RMQ>
void testLargeRmq(std::vector<Elem> values) {
    RMQ rmq(values);
    NlognRMQ<> reference(values);
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Index> idxDist(0, Index(values.size()));
    Index numIterations = Index(std::sqrt(std::sqrt(values.size() + 10'00)));
    for (Index i = 0; i < numIterations; ++i) {
        std::pair<Index, Index> range = std::minmax(idxDist(engine), idxDist(engine));
        if (range.first == range.second || range.first == values.size()) continue;
        testQuery(rmq, reference, range.first, range.second);
    }
    for (Index i = 0; i < numIterations / 10; ++i) {
        Index l = idxDist(engine);
        Index u = idxDist(engine) % 7 + l + 1;
        if (u >= values.size()) { continue; }
        testQuery(rmq, reference, l, u);
    }
}


TEST(LinearSpaceRMQ, BlockBorders) {
    std::vector<Elem> values(randomValues(1 << 18));
    LinearSpaceRMQ<> rmq(values);
    Index blockSize = LinearSpaceRMQ<>::blockSize;
    NlognRMQ<> reference(values);
    for (Index i = 0; i < values.size(); i += blockSize) {
        for (Index j = i + blockSize; j < values.size(); j += blockSize) {
            testQuery(rmq, reference, i, j);
        }
    }
}

TEST(LinearSpaceRMQ, SubBlockBorders) {
    std::vector<Elem> values(randomValues(1 << 16));
    LinearSpaceRMQ<> rmq(values);
    Index subBlockSize = LinearSpaceRMQ<>::subBlockSize;
    NlognRMQ<> reference(values);
    for (Index i = 0, iteration = 0; i < values.size(); i += subBlockSize * ++iteration) {
        for (Index j = i + subBlockSize; j < values.size(); j += subBlockSize) {
            testQuery(rmq, reference, i, j);
        }
    }
}

TEST(LinearSpaceRMQ, AlmostSubBlockBorders) {
    std::vector<Elem> values(randomValues(1 << 16));
    LinearSpaceRMQ<> rmq(values);
    Index subBlockSize = LinearSpaceRMQ<>::subBlockSize;
    NlognRMQ<> reference(values);
    for (Index i = 0, iteration = 0; i < values.size(); i += subBlockSize * ++iteration) {
        for (Index j = i + subBlockSize; j < values.size(); j += subBlockSize) {
            Index l = std::max(int(i + (iteration + j) % 3 - 1), 0);
            Index u = std::min(Index(values.size()), j + (iteration + j) % 3 - 1);
            testQuery(rmq, reference, l, u);
        }
    }
}

TYPED_TEST(AllRmqsTest, Empty) {
    TypeParam rmq = TypeParam();
    ASSERT_EQ(rmq.size(), 0);
}

TYPED_TEST(AllRmqsTest, VerySmall) {
    testSmallRmq<TypeParam>({7});
    testSmallRmq<TypeParam>({3, 6});
    testSmallRmq<TypeParam>({5, 5});
    testSmallRmq<TypeParam>({1, 2, 3});
    testSmallRmq<TypeParam>({1, 2, 1});
    testSmallRmq<TypeParam>({9, 2, 42, 3});
    testSmallRmq<TypeParam>({1, 2, 1, 3});
    testSmallRmq<TypeParam>({0, 2, 3, 3, 100});
}

TYPED_TEST(AllRmqsTest, Small) {
    testSmallRmq<TypeParam>({12, 0, 4, 2, 100, 101});
}


TYPED_TEST(AllRmqsTest, SmallAscending) {
    testSmallRmq<TypeParam>({2, 3, 4, 5, 6, 7, 100, 101, 102, 1000, 1234, 9998, 9999});
}

TYPED_TEST(AllRmqsTest, SmallDescending) {
    testSmallRmq<TypeParam>({9999, 9998, 4567, 3456, 2345, 1234, 1000, 999, 989, 5, 3, 1});
}

TYPED_TEST(AllRmqsTest, SmallRepeating) {
    testSmallRmq<TypeParam>(std::vector<Elem>((1 << 9) + 3, 42));
}

TYPED_TEST(UsefulRmqsTest, PowerOfTwo) {
    testLargeRmq<TypeParam>(randomValues(1 << 16));
}

TYPED_TEST(UsefulRmqsTest, PowerOfTwoMinus1) {
    testLargeRmq<TypeParam>(randomValues((1 << 16) - 1));
}

TYPED_TEST(UsefulRmqsTest, PowerOfTwoPlus1) {
    testLargeRmq<TypeParam>(randomValues((1 << 16) + 1));
}

TYPED_TEST(UsefulRmqsTest, LargeAscending) {
    std::vector<Elem> vec(1 << 18);
    std::iota(vec.begin(), vec.end(), 5);
    testLargeRmq<TypeParam>(std::move(vec));
}

TYPED_TEST(UsefulRmqsTest, LargeDescending) {
    std::vector<Elem> vec(1 << 18);
    std::iota(vec.rbegin(), vec.rend(), 7);
    testLargeRmq<TypeParam>(vec);
}

TYPED_TEST(UsefulRmqsTest, LargeAlmostAscending) {
    std::vector<Elem> vec(1 << 18);
    std::iota(vec.rbegin(), vec.rend(), 7);
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Index> idxDist(0, vec.size() - 1);
    for (Index i = 0; i < 100; ++i) {
        std::swap(vec[idxDist(engine)], vec[idxDist(engine)]);
    }
    testLargeRmq<TypeParam>(vec);
}

TYPED_TEST(UsefulRmqsTest, LargeRepeating) {
    testLargeRmq<TypeParam>(std::vector<Elem>(1 << 20, 54));
}

TYPED_TEST(AllRmqsTest, SmallishRandom) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist(0, 1 << 10);
    testLargeRmq<TypeParam>(randomValues(dist(engine)));
}

TYPED_TEST(UsefulRmqsTest, LargeRandom) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist(0, 1 << 18);
    testLargeRmq<TypeParam>(randomValues(dist(engine)));
}


//template<typename RMQ>
//void testRandomRmq() {
//    auto engine = createRandomEngine();
//    std::uniform_int_distribution<Elem> dist;
//    std::vector<Elem> values(dist(engine) % (1 << 10));
//    for (auto& val: values) {
//        val = dist(engine);
//    }
//    RMQ rmq(values);
//    SimpleRMQ<> reference(values);
//    std::uniform_int_distribution<Index> idxDist(0, Index(values.size()));
//    for (Index i = 0; i < 100'000; ++i) {
//        std::pair<Index, Index> range = std::minmax(idxDist(engine), idxDist(engine));
//        if (range.first == range.second || range.first == values.size()) continue;
//        ASSERT_EQ(rmq(range.first, range.second), reference(range.first, range.second)) << range.first << " " << range.second;
//    }
//}
//TEST(NaiveRmq, Small) {
//    testSmallRmq<NaiveRMQ<>>();
//}
//
//TEST(NaiveRmq, Random) {
//    testRandomRmq<NaiveRMQ<>>();
//}
//
//TEST(NlognRMQ, Small) {
//    testSmallRmq<NlognRMQ<>>();
//}
//
//TEST(NlognRMQ, Random) {
//    testRandomRmq<NlognRMQ<>>();
//}
//
//TEST(LinearSpaceRMQ, Small) {
//    testSmallRmq<LinearSpaceRMQ<>>();
//}

//TEST(CartesianTreeRMQ, Small) {
//    testSmallRmq<CartesianTreeRMQ<>>();
//}
