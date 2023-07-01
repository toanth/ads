
#include "gtest/gtest.h"

#include "../include/common.hpp"

using namespace ads;

TEST(Common, BytesNeededForIndexing) {
    ASSERT_EQ(bytesNeededForIndexing(1), 1);
    ASSERT_EQ(bytesNeededForIndexing(255), 1);
    ASSERT_EQ(bytesNeededForIndexing(256), 1);
    ASSERT_EQ(bytesNeededForIndexing(257), 2);
    ASSERT_EQ(bytesNeededForIndexing(65535), 2);
    ASSERT_EQ(bytesNeededForIndexing(65536), 2);
    ASSERT_EQ(bytesNeededForIndexing(Elem(1) << 32), 4);
    ASSERT_EQ(bytesNeededForIndexing((Elem(1) << 32) + 1), 5);
}

TEST(Common, RoundUpDiv) {
    ASSERT_EQ(roundUpDiv(3, 2), 2);
    ASSERT_EQ(roundUpDiv(0, 4), 0);
    ASSERT_EQ(roundUpDiv(15, 3), 5);
    ASSERT_EQ(roundUpDiv(14, 3), 5);
    ASSERT_EQ(roundUpDiv(16, 3), 6);
}

#ifdef ADS_HAS_CPP20

TEST(Common, Constexpr) {
    static_assert(bytesNeededForIndexing(123) == 1);
    static_assert(bytesNeededForIndexing(257) == 2);
    static_assert(roundUpDiv(53, 8) == 7);
    static_assert(ads::abs(789) == 789);
    static_assert(ads::abs(-1024) == 1024);
}

#endif
