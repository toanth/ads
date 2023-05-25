#include "../include/bitvector.hpp"
#include "gtest/gtest.h"

using namespace ads;

using TestLayout = CacheEfficientLayout;
//using TestLayout = SimpleLayout<>;

TEST(BitvectorConstruction, Sizes) {
    {
        Bitvector<TestLayout> bv(64);
        ASSERT_EQ(bv.numSuperblocks(), 1);
        ASSERT_EQ(bv.sizeInBits(), 64);
        ASSERT_EQ(bv.sizeInElems(), 1);
        for (Index i = 0; i < 100'000; ++i) {
            Bitvector<TestLayout> bv(i);
            ASSERT_EQ(bv.sizeInBits(), i);
            ASSERT_GE(bv.sizeInElems(), (i + 63) / 64);// TODO: Make exact?
            ASSERT_EQ(bv.numSuperblocks(), (bv.sizeInElems() + bv.superblockSize() - 1) / bv.superblockSize());
        }
    }

    Bitvector<TestLayout> bv("1");
    ASSERT_EQ(bv.numSuperblocks(), 1);
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.sizeInElems(), 1);

    std::string text;
    text.reserve(10'000);
    for (Index i = 0; i < 10'000; ++i) {
        Bitvector<TestLayout> bv(text);
        ASSERT_EQ(bv.sizeInBits(), i);
        ASSERT_GE(bv.sizeInElems(), (i + 63) / 64);
        ASSERT_EQ(bv.numSuperblocks(), (bv.sizeInElems() + bv.superblockSize() - 1) / bv.superblockSize());
        text.push_back('1');
    }
}


TEST(BitvectorConstruction, FromStringview) {
    Bitvector<TestLayout> bv("01");
    ASSERT_EQ(bv.sizeInBits(), 2);
    ASSERT_EQ(bv.bit(0), 0);
    ASSERT_EQ(bv.bit(1), 1);
    for (Index i = 2; i < 64; ++i) {
        ASSERT_EQ(bv.bit(i), 0);
    }
    std::string s("101110011101");
    bv = Bitvector<TestLayout>(s);
    ASSERT_EQ(bv.sizeInElems(), 1);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.bit(i), s[i] == '1');
    }
    s = std::string(200, '1');
    s[63] = s[64 + 63] = s[64 * 2 + 63] = '0';
    bv = Bitvector<TestLayout>(s);
    ASSERT_EQ(bv.sizeInElems(), 4);
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.element(i), Elem(-2));
    }
    for (Index i = 0; i < 64; ++i) {
        ASSERT_EQ(bv.bit(64 * 3 + i), i < 200 - 64 * 3);
    }
    bv = Bitvector<TestLayout>(std::string(11, 'e'), 16);
    for (Index i = 0; i < 11 * 4; ++i) {
        ASSERT_EQ(bv.bit(i), i % 4 != 3) << i;
    }
}

TEST(BitvecRank, Only1s) {
    Bitvector<TestLayout> bv("111");
    ASSERT_EQ(bv.element(0), Elem(0b111) << 61);
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
        ASSERT_EQ(bv.rankZero(i), 0) << i;
    }
}

TEST(BitvecRank, Small) {
    Bitvector<TestLayout> bv("110011001100");
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.rankZero(4 * i), 2 * i);
        ASSERT_EQ(bv.rankZero(4 * i + 1), 2 * i);
        ASSERT_EQ(bv.rankZero(4 * i + 2), 2 * i);
        ASSERT_EQ(bv.rankZero(4 * i + 3), 2 * i + 1);
        ASSERT_EQ(bv.rankOne(4 * i), 2 * i);
        ASSERT_EQ(bv.rankOne(4 * i + 1), 2 * i + 1);
        ASSERT_EQ(bv.rankOne(4 * i + 2), 2 * i + 2);
        ASSERT_EQ(bv.rankOne(4 * i + 3), 2 * i + 2);
    }
}

TEST(BitvecRank, OneSuperblock) {
    Bitvector<TestLayout> bv(0);
    std::string s(bv.superblockSize(), '1');
    for (Index i = 0; i < s.size(); i += 2) {
        s[i] = '0';
    }
    bv = Bitvector<TestLayout>(s);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 1) / 2);
        ASSERT_EQ(bv.rankOne(i), i / 2);
    }
}

TEST(BitvecRank, ManySuperblocks) {
    Bitvector<TestLayout> bv(0);
    std::string s(bv.superblockSize() * 7 + 123456, '1');
    for (Index i = 0; i < s.size(); i += 3) {
        s[i] = '0';
    }
    bv = Bitvector<TestLayout>(s);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 2) / 3);
        ASSERT_EQ(bv.rankOne(i), (2 * i) / 3);
    }
}

TEST(BitvecSelect, Only1s) {
    Bitvector<TestLayout> bv("1111");
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.selectOne(i), i) << i;
        ASSERT_EQ(bv.selectZero(i), -1) << i;
    }
}

TEST(BitvecSelect, Small) {
    std::string s("1100000000000000000000000000001");
    Bitvector<TestLayout> bv(s);
    ASSERT_EQ(bv.selectOne(0), 0);
    ASSERT_EQ(bv.selectOne(1), 1);
    ASSERT_EQ(bv.selectOne(2), s.size() - 1);
    for (Index i = 0; i < s.size() - 3; ++i) {
        ASSERT_EQ(bv.selectZero(i), i + 2);
    }
}

TEST(BitvecSelect, Large) {
    std::string s(15, 'c');
    Bitvector<TestLayout> bv(s, 16);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.selectZero(2 * i), 4 * i + 2);
        ASSERT_EQ(bv.selectZero(2 * i + 1), 4 * i + 3);
        ASSERT_EQ(bv.selectOne(2 * i), 4 * i);
        ASSERT_EQ(bv.selectOne(2 * i + 1), 4 * i + 1);
    }
}


TEST(Bitvector, EmptyOrOneElem) {
    Bitvector<TestLayout> bv(0);
    ASSERT_EQ(bv.sizeInBits(), 0);
    ASSERT_EQ(bv.sizeInElems(), 0);
    bv = Bitvector<TestLayout>("1");
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.sizeInElems(), 1);
    ASSERT_EQ(bv.rankOne(0), 0);
    ASSERT_EQ(bv.rankZero(0), 0);
}
