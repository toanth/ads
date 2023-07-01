
#include "gtest/gtest.h"

#include "../include/bit.hpp"

using namespace ads;


TEST(Bit, IntLog2) {
    ASSERT_EQ(intLog2(1u), 0);
    ASSERT_EQ(intLog2(std::uint64_t(2)), 1);
    ASSERT_EQ(intLog2(std::uint32_t(3)), 1);
    ASSERT_EQ(intLog2(std::uint8_t(4)), 2);
    ASSERT_EQ(intLog2(std::uint8_t(5)), 2);
    ASSERT_EQ(intLog2(std::uint16_t(123)), 6);
    ASSERT_EQ(intLog2(std::uint8_t(255)), 7);
    ASSERT_EQ(intLog2(0x19ull), 4);
}

TEST(Bit, RoundUpLog2) {
    ASSERT_EQ(roundUpLog2(1ull), 0);
    ASSERT_EQ(roundUpLog2(2ul), 1);
    ASSERT_EQ(roundUpLog2(3u), 2);
    ASSERT_EQ(roundUpLog2(4u), 2);
    ASSERT_EQ(roundUpLog2(5u), 3);
    ASSERT_EQ(roundUpLog2(std::uint8_t(123)), 7);
    ASSERT_EQ(roundUpLog2(std::uint8_t(255)), 8);
}


TEST(Bit, Popcount) {
    ASSERT_EQ(popcount(0u), 0);
    ASSERT_EQ(popcount(std::uint16_t(1)), 1);
    ASSERT_EQ(popcount(std::uint8_t(2)), 1);
    ASSERT_EQ(popcount(3ul), 2);
    ASSERT_EQ(popcount(std::uint8_t(255)), 8);
    ASSERT_EQ(popcount(std::uint16_t(256)), 1);
    ASSERT_EQ(popcount((1ull << 16) - 12), 13);
    ASSERT_EQ(popcount(Elem(-1)), 64);
    ASSERT_EQ(popcount(Elem(-2)), 63);
}

TEST(Bit, PopcountFallback) {
    ASSERT_EQ(popcountFallback(0u), 0);
    ASSERT_EQ(popcountFallback(std::uint16_t(1)), 1);
    ASSERT_EQ(popcountFallback(std::uint8_t(2)), 1);
    ASSERT_EQ(popcountFallback(3ul), 2);
    ASSERT_EQ(popcountFallback(std::uint8_t(255)), 8);
    ASSERT_EQ(popcountFallback(std::uint16_t(256)), 1);
    ASSERT_EQ(popcountFallback((1ull << 16) - 12), 13);
    ASSERT_EQ(popcountFallback(Elem(-1)), 64);
    ASSERT_EQ(popcountFallback(Elem(-2)), 63);
}

TEST(Bit, ReverseBits) {
    ASSERT_EQ(reverseBits(std::uint8_t(0)), 0);
    ASSERT_EQ(reverseBits(std::uint8_t(1)), 128);
    ASSERT_EQ(reverseBits(std::uint64_t(1)), std::uint64_t(1) << 63);
    ASSERT_EQ(reverseBits(std::uint16_t(3)), (1 << 15) + (1 << 14));
    ASSERT_EQ(reverseBits(std::uint64_t(-2)), std::uint64_t(-1) / 2);
    ASSERT_EQ(reverseBits(0x1248'ffff'0110'fedcull), 0x3b7f'0880'ffff'1248ull);
    std::uint64_t v = 18374161034817641865ull;
    ASSERT_EQ(reverseBits(reverseBits(v)), v);
}

TEST(Bit, CountTrailingZeros) {
    ASSERT_EQ(countTrailingZeros(1), 0);
    ASSERT_EQ(countTrailingZeros(2), 1);
    ASSERT_EQ(countTrailingZeros(3), 0);
    ASSERT_EQ(countTrailingZeros(Elem(1) << 32), 32);
    ASSERT_EQ(countTrailingZeros(Elem(-2)), 1);
    ASSERT_EQ(countTrailingZeros(Elem(-1)), 0);
}

TEST(Bit, ElemSelect) {
    ASSERT_EQ(elemSelect(1, 0), 0);
    ASSERT_EQ(elemSelect(2, 0), 1);
    ASSERT_EQ(elemSelect(3, 0), 0);
    ASSERT_EQ(elemSelect(3, 1), 1);
    ASSERT_EQ(elemSelect((Elem(1) << 15) + 3, 2), 15);
    ASSERT_EQ(elemSelect(Elem(-1), 7), 7);
    ASSERT_EQ(elemSelect(Elem(-1), 63), 63);
}

TEST(Bit, ConstevalElemSelect) {
    static_assert(constevalElemSelect(1, 0) == 0);
    static_assert(constevalElemSelect(2, 0) == 1);
    static_assert(constevalElemSelect(3, 0) == 0);
    static_assert(constevalElemSelect(3, 1) == 1);
    static_assert(constevalElemSelect((Elem(1) << 15) + 3, 2) == 15);
    static_assert(constevalElemSelect(Elem(-1), 7) == 7);
    static_assert(constevalElemSelect(Elem(-1), 63) == 63);
}

TEST(Bit, PrecomputeBitSelectTable4Bits) {
    auto table = precomputeBitSelectTable<4>();
    ASSERT_EQ(table.size(), 16);
    ASSERT_EQ(table[0].size(), 4);
    ASSERT_EQ(table[2][0], 1);
    ASSERT_EQ(table[15][3], 3);
    ASSERT_EQ(table[11][1], 1);
    ASSERT_EQ(table[13][0], 0);
}

TEST(Bit, PrecomputeBitSelectTable) {
    auto table = precomputeBitSelectTable();
    ASSERT_EQ(table.size(), 256);
    ASSERT_EQ(table[0].size(), 8);
    ASSERT_EQ(table[1][0], 0);
    ASSERT_EQ(table[2][0], 1);
    ASSERT_EQ(table[3][0], 0);
    ASSERT_EQ(table[3][1], 1);
    ASSERT_EQ(table[11][2], 3);
    ASSERT_EQ(table[(1 << 7) + 4][1], 7);
    ASSERT_EQ(table[(1 << 6) + 9][2], 6);
    for (Index i = 0; i < 8; ++i) {
        ASSERT_EQ(table[255][i], i);
    }
}

TEST(Bit, ElemSelectWithTable) {
    ASSERT_EQ(elemSelectWithTable(1, 0), 0);
    ASSERT_EQ(elemSelectWithTable(2, 0), 1);
    ASSERT_EQ(elemSelectWithTable(3, 0), 0);
    ASSERT_EQ(elemSelectWithTable(3, 1), 1);
    ASSERT_EQ(elemSelectWithTable((Elem(1) << 15) + 3, 2), 15);
    ASSERT_EQ(elemSelectWithTable(Elem(-1), 7), 7);
    ASSERT_EQ(elemSelectWithTable(Elem(-1), 63), 63);
}


#ifdef ADS_HAS_CPP20

TEST(Bit, Constexpr) {
    static_assert(intLog2(987u) == 9);
    static_assert(roundUpLog2(42ul) == 6);
    static_assert(popcount(0xabcu) == 2 + 3 + 2);
    static_assert(reverseBits(0x0123u) == (0xc480u << 16));
    static_assert(countTrailingZeros((unsigned char)14) == 1);
    static_assert(elemSelect(0x1234ull, 3) == 9);
    static_assert(constevalElemSelect(1243619ull, 12) == elemSelect(1243619, 12));
    static_assert(elemSelectWithTable(127ull, 4) == 4);
}

#endif
