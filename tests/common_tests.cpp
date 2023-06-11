
#include "gtest/gtest.h"

#include "../include/common.hpp"

using namespace ads;

TEST(bytesNeededForIndexing, General) {
    ASSERT_EQ(bytesNeededForIndexing(1), 1);
    ASSERT_EQ(bytesNeededForIndexing(255), 1);
    ASSERT_EQ(bytesNeededForIndexing(256), 1);
    ASSERT_EQ(bytesNeededForIndexing(257), 2);
    ASSERT_EQ(bytesNeededForIndexing(65535), 2);
    ASSERT_EQ(bytesNeededForIndexing(65536), 2);
    ASSERT_EQ(bytesNeededForIndexing(Elem(1) << 32), 4);
    ASSERT_EQ(bytesNeededForIndexing((Elem(1) << 32) + 1), 5);
}

TEST(Log2, General) {
    ASSERT_EQ(log2(1u), 0);
    ASSERT_EQ(log2(std::uint64_t(2)), 1);
    ASSERT_EQ(log2(std::uint32_t(3)), 1);
    ASSERT_EQ(log2(std::uint8_t(4)), 2);
    ASSERT_EQ(log2(std::uint8_t(5)), 2);
    ASSERT_EQ(log2(std::uint16_t(123)), 6);
    ASSERT_EQ(log2(std::uint8_t(255)), 7);
    ASSERT_EQ(log2(0x19ull), 4);
}

TEST(RoundUpLog2, General) {
    ASSERT_EQ(roundUpLog2(1ull), 0);
    ASSERT_EQ(roundUpLog2(2ul), 1);
    ASSERT_EQ(roundUpLog2(3u), 2);
    ASSERT_EQ(roundUpLog2(4u), 2);
    ASSERT_EQ(roundUpLog2(5u), 3);
    ASSERT_EQ(roundUpLog2(std::uint8_t(123)), 7);
    ASSERT_EQ(roundUpLog2(std::uint8_t(255)), 8);
}


TEST(Popcount, General) {
    ASSERT_EQ(popcount(0), 0);
    ASSERT_EQ(popcount(std::uint16_t(1)), 1);
    ASSERT_EQ(popcount(std::uint8_t(2)), 1);
    ASSERT_EQ(popcount(3), 2);
    ASSERT_EQ(popcount(255), 8);
    ASSERT_EQ(popcount(256), 1);
    ASSERT_EQ(popcount((1 << 16) - 12), 13);
    ASSERT_EQ(popcount(Elem(-1)), 64);
    ASSERT_EQ(popcount(Elem(-2)), 63);
}

TEST(ReverseBits, General) {
    ASSERT_EQ(reverseBits(std::uint8_t(0)), 0);
    ASSERT_EQ(reverseBits(std::uint8_t(1)), 128);
    ASSERT_EQ(reverseBits(std::uint64_t(1)), std::uint64_t(1) << 63);
    ASSERT_EQ(reverseBits(std::uint16_t(3)), (1 << 15) + (1 << 14));
    ASSERT_EQ(reverseBits(std::uint64_t(-2)), std::uint64_t(-1) / 2);
    ASSERT_EQ(reverseBits(0x1248'ffff'0110'fedcull), 0x3b7f'0880'ffff'1248ull);
    std::uint64_t v = 18374161034817641865ull;
    ASSERT_EQ(reverseBits(reverseBits(v)), v);
}

TEST(RoundUpDiv, General) {
    ASSERT_EQ(roundUpDiv(3, 2), 2);
    ASSERT_EQ(roundUpDiv(0, 4), 0);
    ASSERT_EQ(roundUpDiv(15, 3), 5);
    ASSERT_EQ(roundUpDiv(14, 3), 5);
    ASSERT_EQ(roundUpDiv(16, 3), 6);
}
