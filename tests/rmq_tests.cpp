
#include "../include/linear_space_rmq.hpp"
#include "../include/naive_rmq.hpp"
#include "../include/nlogn_rmq.hpp"
#include "../include/succinct_rmq.hpp"
#include "gtest/gtest.h"

using namespace ads;


struct NoCheck {
    template<typename... Args>
    constexpr bool operator()(Args&&...) noexcept {
        return true;
    }
};


template<ADS_RMQ_CONCEPT RmqType>
class AllRmqsTest : public ::testing::Test {};

template<ADS_RMQ_CONCEPT RmqType>
class UsefulRmqsTest : public ::testing::Test {};

using AllRmqTypes = ::testing::Types<NaiveRMQ<>, NLogNRmq<>, LinearSpaceRMQ<>, SuccinctRMQ<>>;
TYPED_TEST_SUITE(AllRmqsTest, AllRmqTypes);

using UsefulRmqTypes = ::testing::Types<NLogNRmq<>, LinearSpaceRMQ<>, SuccinctRMQ<>>;
TYPED_TEST_SUITE(UsefulRmqsTest, UsefulRmqTypes);


void fillWithRandomValues(std::vector<Elem>& v, Elem maxVal) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist(maxVal);
    for (Elem& value : v) {
        value = dist(engine);
    }
}

std::vector<Elem> randomValues(Index size, Elem maxVal = std::numeric_limits<Elem>::max()) {
    std::vector<Elem> res(size);
    fillWithRandomValues(res, maxVal);
    return res;
}

template<ADS_RMQ_CONCEPT TestRmq, ADS_RMQ_CONCEPT ReferenceRmq, typename CheckAnswer = NoCheck>
void testQuery(const TestRmq& testRmq, const ReferenceRmq& referenceRmq, Index lower, Index upper,
        CheckAnswer check = CheckAnswer()) {
    Index testeeAnswer = testRmq(lower, upper);
    Index referenceAnswer = referenceRmq(lower, upper);
    ASSERT_EQ(referenceRmq[testeeAnswer], referenceRmq[referenceAnswer])
            << TestRmq::name << ", answers (test/ref) " << testeeAnswer << " " << referenceAnswer << ", interval "
            << lower << " " << upper << ", size " << referenceRmq.size();
    ASSERT_TRUE(check(lower, upper, testeeAnswer));
}

template<typename RMQ>
void testSmallRmq(std::vector<Elem> values) {
    RMQ rmq(values);
    SimpleRMQ<> reference(values);
    for (Index i = 0; i < ssize(values); ++i) {
        for (Index j = i + 1; j <= ssize(values); ++j) {
            testQuery(rmq, reference, i, j);
        }
    }
}

template<typename RMQ, typename CheckAnswer = NoCheck>
void testLargeRmq(std::vector<Elem> values, CheckAnswer check = CheckAnswer()) {
    RMQ rmq(values);
    NLogNRmq<> reference(values);
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Index> idxDist(0, ssize(values));
    Index numIterations = ssize(values) / 10 + 10 * Index(std::sqrt(ssize(values))) + 10;
    for (Index i = 0; i < numIterations; ++i) {
        std::pair<Index, Index> range = std::minmax(idxDist(engine), idxDist(engine));
        if (range.first == range.second) continue;
        testQuery(rmq, reference, range.first, range.second);
        if (testing::Test::HasFatalFailure()) {
            return;
        }
    }
    for (Index i = 0; i < numIterations / 4; ++i) {
        Index l = idxDist(engine);
        Index u = idxDist(engine) % 7 + l + 1;
        if (u >= ssize(values)) {
            continue;
        }
        testQuery(rmq, reference, l, u, check);
        if (testing::Test::HasFatalFailure()) {
            return;
        }
    }
}


TEST(LinearSpaceRMQ, BlockBorders) {
    std::vector<Elem> values(randomValues(1 << 18));
    LinearSpaceRMQ<> rmq(values);
    Index blockSize = LinearSpaceRMQ<>::blockSize;
    NLogNRmq<> reference(values);
    for (Index i = 0; i < ssize(values); i += blockSize) {
        for (Index j = i + blockSize; j < ssize(values); j += blockSize) {
            testQuery(rmq, reference, i, j);
        }
    }
}

TEST(LinearSpaceRMQ, SubBlockBorders) {
    std::vector<Elem> values(randomValues(1 << 16));
    LinearSpaceRMQ<> rmq(values);
    Index subBlockSize = LinearSpaceRMQ<>::subBlockSize;
    NLogNRmq<> reference(values);
    for (Index i = 0, iteration = 0; i < ssize(values); i += subBlockSize * ++iteration) {
        for (Index j = i + subBlockSize; j < ssize(values); j += subBlockSize) {
            testQuery(rmq, reference, i, j);
        }
    }
}

TEST(LinearSpaceRMQ, AlmostSubBlockBorders) {
    std::vector<Elem> values(randomValues(1 << 16));
    LinearSpaceRMQ<> rmq(values);
    Index subBlockSize = LinearSpaceRMQ<>::subBlockSize;
    NLogNRmq<> reference(values);
    for (Index i = 0, iteration = 0; i < ssize(values); i += subBlockSize * ++iteration) {
        for (Index j = i + subBlockSize; j < ssize(values); j += subBlockSize) {
            Index l = std::max(int(i + (iteration + j) % 3 - 1), 0);
            Index u = std::min(ssize(values), j + (iteration + j) % 3 - 1);
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
    testSmallRmq<TypeParam>({20, 8, 1, 15, 0});
}

TYPED_TEST(AllRmqsTest, Small) {
    testSmallRmq<TypeParam>({12, 0, 4, 2, 100, 101, 79, 5, 2, 100});
}

TYPED_TEST(AllRmqsTest, SuccinctRmqPaperExample) {
    TypeParam rmq{5, 3, 4, 3, 4, 5, 1, 3, 2, 4, 2, 5, 3, 5, 5, 4};
    ASSERT_EQ(rmq(6, 9), 6);
    ASSERT_EQ(rmq(9, 15), 10);
}

TYPED_TEST(AllRmqsTest, SmallAscending) {
    testSmallRmq<TypeParam>({2, 3, 4, 5, 6, 7, 100, 101, 102, 1000, 1234, 9998, 9999});
}

TYPED_TEST(AllRmqsTest, SmallDescending) {
    testSmallRmq<TypeParam>({9999, 9998, 4567, 3456, 2345, 1234, 1000, 999, 989, 5, 3, 1});
}

TYPED_TEST(AllRmqsTest, SmallRepeating) {
    testSmallRmq<TypeParam>(std::vector<Elem>((1 << 8) + 13, 42));
}

TYPED_TEST(AllRmqsTest, SuccinctBlockBorder) {
    std::vector<Elem> vec(260);
    for (Index i = 0; i < ssize(vec); ++i) {
        vec[i] = i;
    }
    TypeParam rmq(vec);
    ASSERT_EQ(rmq(240, 260), 240);
    for (Elem& v : vec) {
        v = vec.back() - v;
    }
    rmq = TypeParam(vec);
    ASSERT_EQ(rmq(240, 260), 259);
    vec.clear();
    for (Index i = 0; i < 240; ++i) {
        vec.push_back(240 - i);
    }
    for (Index i = 240; i < 350; ++i) {
        vec.push_back(1000 - i);
    }
    rmq = TypeParam(vec);
    ASSERT_EQ(rmq(240, 260), 259);
    ASSERT_EQ(rmq(239, 260), 239);
    ASSERT_EQ(rmq(200, 350), 239);
    ASSERT_EQ(rmq(240, 241), 240);
    ASSERT_EQ(rmq(239, 240), 239);
}

TYPED_TEST(AllRmqsTest, 1000Elements) {
    std::vector<Elem> vec(1000);
    for (Index i = 0; i * 2 < ssize(vec); ++i) {
        vec[2 * i] = i + 1'000'200;
        vec[2 * i + 1] = (i % 2 == 0 ? 1 : -1) * i * i + 1'000'000;
    }
    testLargeRmq<TypeParam>(vec);
}

TYPED_TEST(AllRmqsTest, SmallManyRandomDuplicates) {
    testSmallRmq<TypeParam>(randomValues(301, 7));
}

TYPED_TEST(AllRmqsTest, SmallishRandom) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Index> dist(0, 1 << 10);
    testLargeRmq<TypeParam>(randomValues(dist(engine)));
}

TYPED_TEST(UsefulRmqsTest, PowerOfTwo) {
    testLargeRmq<TypeParam>(randomValues(1 << 16));
}

TYPED_TEST(UsefulRmqsTest, PowerOfTwoMinus1) {
    testLargeRmq<TypeParam>(randomValues((1 << 16) - 1));
}

// This test exists because the bitvector of the succinct rmq stores 2 (n+1) bits
TYPED_TEST(UsefulRmqsTest, PowerOfTwoMinus2) {
    testLargeRmq<TypeParam>(randomValues((1 << 16) - 2));
}

TYPED_TEST(UsefulRmqsTest, PowerOfTwoPlus1) {
    testLargeRmq<TypeParam>(randomValues((1 << 16) + 1));
}

TYPED_TEST(UsefulRmqsTest, LargeAscending) {
    std::vector<Elem> vec((1 << 18) - 7);
    std::iota(vec.begin(), vec.end(), 5);
    auto takeLowest = [](Index low, Index, Index answer) { return answer == low; };
    testLargeRmq<TypeParam>(std::move(vec), takeLowest);
}

TYPED_TEST(UsefulRmqsTest, LargeDescending) {
    std::vector<Elem> vec((1 << 18) + 23);
    std::iota(vec.rbegin(), vec.rend(), 7);
    auto takeHighest = [](Index, Index high, Index answer) { return answer == high - 1; };
    testLargeRmq<TypeParam>(vec, takeHighest);
}

TYPED_TEST(UsefulRmqsTest, LargeAlmostAscending) {
    std::vector<Elem> vec((1 << 18) - 31);
    std::iota(vec.rbegin(), vec.rend(), 7);
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Index> idxDist(0, ssize(vec) - 1);
    for (Index i = 0; i < 100; ++i) {
        std::swap(vec[idxDist(engine)], vec[idxDist(engine)]);
    }
    testLargeRmq<TypeParam>(vec);
}

TYPED_TEST(UsefulRmqsTest, LargeRepeating) {
    testLargeRmq<TypeParam>(std::vector<Elem>((1 << 17) - 4, 54));
}

TYPED_TEST(UsefulRmqsTest, LargeRandom) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist(0, 1 << 18);
    testLargeRmq<TypeParam>(randomValues(Index(dist(engine)), 100));
}


#ifdef ADS_HAS_CPP20
TYPED_TEST(AllRmqsTest, Constexpr) {
    static_assert(TypeParam{}.size() == 0);
    static_assert(TypeParam{3}(0, 1) == 0);
    static_assert(TypeParam{0, 23, 2, 40, 4, 2, 1}(2, 5) == 2);
    static constexpr Index size = std::same_as<TypeParam, NaiveRMQ<>>     ? 50 :
                                  std::same_as<TypeParam, NLogNRmq<Elem>> ? 500 :
                                                                            1234;
    constexpr auto testVec = []() {
        std::vector<Elem> res(size);
        std::iota(res.begin(), res.end(), 987654);
        res[42] = 123456;
        return res;
    };
    static_assert(TypeParam(testVec())(12, size - 1) == 42);
}
#endif
