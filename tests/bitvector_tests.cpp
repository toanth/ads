#include "../include/bitvector.hpp"
#include "gtest/gtest.h"

using namespace ads;

using SmallSuperBlocks = SimpleLayout<1, (1 << 8) / 64, std::uint32_t>;
using NonsensicalSuperBlocks = SimpleLayout<2, (1 << 14) / 64, std::uint64_t>;
using OneBlockPerSuperBlock = SimpleLayout<64, 64>;
using LargeBlockCounts = SimpleLayout<4, (1 << 16)>;

template<ADS_LAYOUT_CONCEPT Layout>
class ManyBitvecLayoutsTest : public ::testing::Test {};

template<ADS_LAYOUT_CONCEPT Layout>
class EfficientBitvecLayoutsTest : public ::testing::Test {};

using ManyLayouts
        = ::testing::Types<CacheEfficientLayout, SimpleLayout<>, SmallSuperBlocks, NonsensicalSuperBlocks, OneBlockPerSuperBlock, LargeBlockCounts>;
using EfficientLayouts = ::testing::Types<CacheEfficientLayout, SimpleLayout<>>;
TYPED_TEST_SUITE(ManyBitvecLayoutsTest, ManyLayouts);
TYPED_TEST_SUITE(EfficientBitvecLayoutsTest, EfficientLayouts);

TYPED_TEST(ManyBitvecLayoutsTest, ConstructionSizes) {
    {
        Allocation<> allocation(Bitvector<TypeParam>::allocatedSizeInElems(1));
        Bitvector<TypeParam> bv(64, allocation.memory());
        ASSERT_EQ(bv.numSuperBlocks(), 1);
        ASSERT_EQ(bv.sizeInBits(), 64);
        ASSERT_EQ(bv.sizeInElems(), 1);
        ASSERT_EQ(bv.allocatedSizeInElems(), allocation.size());
        for (Index i = 0; i < 100'000; ++i) {
            bv = Bitvector<TypeParam>(i);
            ASSERT_EQ(bv.sizeInBits(), i);
            ASSERT_EQ(bv.bitView().size(), bv.sizeInBits());
            ASSERT_GE(bv.sizeInElems(), (i + 63) / 64); // TODO: Make exact?
            ASSERT_EQ(bv.elemView().size(), bv.sizeInElems());
            ASSERT_EQ(bv.numSuperBlocks(), (bv.sizeInElems() + bv.numElemsInSuperBlock() - 1) / bv.numElemsInSuperBlock());
        }
    }

    Bitvector<TypeParam> bv("1");
    ASSERT_EQ(bv.numSuperBlocks(), 1);
    ASSERT_EQ(bv.numBlocks(), 1);
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.sizeInElems(), 1);

    std::string text;
    text.reserve(10'000);
    for (Index i = 0; i < 10'000; ++i) {
        Bitvector<TypeParam> bv(text);
        ASSERT_EQ(bv.sizeInBits(), i);
        ASSERT_GE(bv.sizeInElems(), (i + 63) / 64);
        ASSERT_EQ(bv.numSuperBlocks(), (bv.sizeInElems() + bv.numElemsInSuperBlock() - 1) / bv.numElemsInSuperBlock());
        text.push_back('1');
    }
}

TYPED_TEST(ManyBitvecLayoutsTest, ConstructionElements) {
    Bitvector<TypeParam> bv(1);
    bv.setElem(0, 1);
    bv.buildRankMetadata(0);
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.getElem(0), 1);
    ASSERT_EQ(bv.getBit(0), 1);
    ASSERT_EQ(*bv.bitView().begin(), 1);
    ASSERT_EQ(bv.toString(), "1");
    ASSERT_EQ(bv, Bitvector<>("1"));
    ASSERT_GT(bv, Bitvector<>("0"));
    ASSERT_EQ(bv.numZeros(), 0);
    ASSERT_EQ(bv.numOnes(), 1);
    bv = Bitvector<TypeParam>(64);
    bv.setElem(0, Elem(1) << 63);
    bv.buildRankMetadata(0);
    ASSERT_EQ(bv.sizeInBits(), 64);
    ASSERT_EQ(bv.getBit(0), 0);
    ASSERT_EQ(bv.getBit(63), 1);
    ASSERT_FALSE(*bv.bitView().begin());
    ASSERT_GT(bv, Bitvector<>("0"));
    ASSERT_EQ(bv.numZeros(), 63);
    ASSERT_EQ(bv.numOnes(), 1);
}

TYPED_TEST(ManyBitvecLayoutsTest, ConstructionFromStringview) {
    Bitvector<TypeParam> bv("01");
    ASSERT_EQ(bv.sizeInBits(), 2);
    ASSERT_EQ(bv.getBit(0), 0);
    ASSERT_EQ(bv.getBit(1), 1);
    for (Index i = 2; i < 64; ++i) {
        ASSERT_EQ(bv.getBit(i), 0);
    }
    ASSERT_EQ(bv.toString(), "01");
    ASSERT_EQ(bv.numZeros(), bv.numOnes());
    std::string s("101110011101");
    bv = Bitvector<TypeParam>(s);
    ASSERT_EQ(bv.sizeInElems(), 1);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.getBit(i), s[i] == '1');
    }
    ASSERT_EQ(bv.numZeros(), 4);
    s = std::string(200, '1');
    s[63] = s[64 + 63] = s[64 * 2 + 63] = '0';
    Allocation<> allocation(Bitvector<TypeParam>::allocatedSizeInElems(roundUpDiv(s.size(), 64)));
    bv = Bitvector<TypeParam>(s, 2, allocation.memory());
    ASSERT_EQ(bv.sizeInElems(), 4);
    ASSERT_EQ(bv.allocatedSizeInElems(), allocation.size());
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.getElem(i), Elem(-1) / 2);
    }
    for (Index i = 0; i < 64; ++i) {
        ASSERT_EQ(bv.getBit(64 * 3 + i), i < 200 - 64 * 3);
    }
    bv = Bitvector<TypeParam>(std::string(11, 'e'), 16);
    for (Index i = 0; i < 11 * 4; ++i) {
        ASSERT_EQ(bv.getBit(i), i % 4 != 3) << i;
    }
}

TYPED_TEST(ManyBitvecLayoutsTest, RankOnly1s) {
    Bitvector<TypeParam> bv("111");
    ASSERT_EQ(bv.getElem(0), Elem(0b111));
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
        ASSERT_EQ(bv.rankZero(i), 0) << i;
    }
    ASSERT_EQ(bv.numZeros(), 0);
}

TYPED_TEST(ManyBitvecLayoutsTest, RankSmall) {
    Bitvector<TypeParam> bv("110011001100");
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

TYPED_TEST(ManyBitvecLayoutsTest, RankOneSuperblock) {
    Bitvector<TypeParam> bv(0);
    std::string s(bv.superBlockSize(), '1');
    for (Index i = 0; i < s.size(); i += 2) {
        s[i] = '0';
    }
    bv = Bitvector<TypeParam>(s);
    ASSERT_EQ(bv.toString(), s);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 1) / 2) << i;
        ASSERT_EQ(bv.rankOne(i), i / 2) << i;
    }
    ASSERT_EQ(bv.numZeros(), bv.sizeInBits() / 2);
}

TYPED_TEST(ManyBitvecLayoutsTest, RankManySuperblocks) {
    Bitvector<TypeParam> bv(0);
    std::string s(bv.superBlockSize() * 7 + 123456, '1');
    for (Index i = 0; i < s.size(); i += 3) {
        s[i] = '0';
    }
    bv = Bitvector<TypeParam>(s);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 2) / 3) << i;
        ASSERT_EQ(bv.rankOne(i), (2 * i) / 3) << i;
    }
}

TYPED_TEST(ManyBitvecLayoutsTest, SelectOnly1s) {
    Bitvector<TypeParam> bv("1111");
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.selectOne(i), i) << i;
        //        ASSERT_EQ(bv.selectZero(i), -1) << i;
    }
    bv = Bitvector<TypeParam>(std::string(200'000, '1'));
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
    }
}

TYPED_TEST(ManyBitvecLayoutsTest, SelectSmall) {
    std::string s("1100000000000000000000000000001");
    Bitvector<TypeParam> bv(s);
    ASSERT_EQ(bv.selectOne(0), 0);
    ASSERT_EQ(bv.selectOne(1), 1);
    ASSERT_EQ(bv.selectOne(2), s.size() - 1);
    for (Index i = 0; i < s.size() - 3; ++i) {
        ASSERT_EQ(bv.selectZero(i), i + 2);
    }
    bv = Bitvector<TypeParam>("0111");
    ASSERT_EQ(bv.selectOne(0), 1);
    ASSERT_EQ(bv.selectOne(1), 2);
    ASSERT_EQ(bv.selectZero(0), 0);
}

TYPED_TEST(ManyBitvecLayoutsTest, SelectLarge) {
    std::string s(12345, 'c');
    Bitvector<TypeParam> bv(s, 16);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.selectZero(2 * i), 4 * i + 2);
        ASSERT_EQ(bv.selectZero(2 * i + 1), 4 * i + 3);
        ASSERT_EQ(bv.selectOne(2 * i), 4 * i);
        ASSERT_EQ(bv.selectOne(2 * i + 1), 4 * i + 1);
    }
}

TYPED_TEST(ManyBitvecLayoutsTest, Select2SuperBlocksPlus2) {
    Bitvector<TypeParam> bv;
    bv = Bitvector<TypeParam>(bv.superBlockSize() * 2 + 2, 1);
    ASSERT_EQ(bv.rankOne(bv.sizeInBits() - 1), bv.sizeInElems());
    for (Index i = 0; i < bv.sizeInBits() - bv.sizeInElems(); ++i) {
        ASSERT_EQ(bv.selectZero(i), i + i / 63 + 1) << i << " " << bv.sizeInBits();
    }
    for (Index i = 0; i < bv.sizeInElems(); ++i) {
        ASSERT_EQ(bv.selectOne(i), i * 64) << i << " " << bv.sizeInBits();
    }
}


TYPED_TEST(ManyBitvecLayoutsTest, EmptyOrOneElem) {
    Bitvector<TypeParam> bv(0);
    ASSERT_EQ(bv.sizeInBits(), 0);
    ASSERT_EQ(bv.sizeInElems(), 0);
    bv = Bitvector<TypeParam>("1");
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.sizeInElems(), 1);
    ASSERT_EQ(bv.rankOne(0), 0);
    ASSERT_EQ(bv.rankZero(0), 0);
}

TYPED_TEST(EfficientBitvecLayoutsTest, PowerOfTwo) {
    for (Index i = 1; i <= (Elem(1) << 26); i *= 4) {
        Bitvector<TypeParam> bv(i, Elem(0));
        for (Index j = 0; j < bv.numElems() - 1; ++j) {
            bv.setElem(j, 0xaaaa'aaaa'aaaa'aaaaull);
        }
        for (Index j = (bv.numElems() - 1) * 64 + 1; j < bv.sizeInBits(); j += 2) {
            bv.setBit(j, true);
        }
        for (Index j = 0; j < bv.numSuperBlocks(); ++j) {
            bv.buildRankMetadata(j);
        }
        for (Index j = 1; j < bv.sizeInBits(); j = j * 3 - 1) {
            ASSERT_EQ(bv.rankZero(j), (j + 1) / 2) << j << " " << i;
            ASSERT_EQ(bv.selectOne(j / 2), j / 2 * 2 + 1) << j << " " << i;
        }
        ASSERT_EQ(bv.numOnes(), bv.sizeInBits() / 2);
    }
}

TYPED_TEST(ManyBitvecLayoutsTest, Random) {
    auto engine = createRandomEngine();
    auto dist = std::uniform_int_distribution<Index>(0, 400'000);
    std::string str(dist(engine), '0');
    for (char& c : str) {
        if (dist(engine) & 1) {
            c = '1';
        }
    }
    std::vector<std::uint32_t> results(str.size());
    results[0] = 0;
    for (Index i = 1; i < str.size(); ++i) {
        results[i] = results[i - 1] + (str[i - 1] == '0');
    }
    Bitvector<TypeParam> bv(str);
    Index numOnes = bv.rankOne(bv.sizeInBits() - 1);
    if (bv.getBit(bv.sizeInBits() - 1)) {
        ++numOnes;
    }
    ASSERT_EQ(numOnes, bv.numOnes());
    ASSERT_EQ(bv.numZeros(), results.back() + (bv.getBit(bv.sizeInBits() - 1) ? 0 : 1));
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.rankZero(i), results[i]) << i << " " << bv.sizeInBits();
        Index rank = bv.rankOne(i);
        if (rank < numOnes) {
            ASSERT_GE(bv.selectOne(rank), i) << bv.sizeInBits();
        } else {
            ASSERT_EQ(rank, numOnes);
        }
        if (i < numOnes) {
            ASSERT_EQ(bv.rankOne(bv.selectOne(i)), i) << bv.sizeInBits();
        }
    }
}

TYPED_TEST(ManyBitvecLayoutsTest, RandomLongRuns) {
    auto engine = createRandomEngine();
    auto dist = std::uniform_real_distribution<double>(-16.0, 16.0);
    std::string str;
    while (str.size() < 300'000) {
        double randomVal = dist(engine);
        Index len = Index(std::log2(std::abs(randomVal)));
        char c = randomVal > 0.0 ? '1' : '0';
        for (Index i = 0; i < len; ++i) {
            str += c;
        }
    }
    std::vector<std::uint32_t> results(str.size());
    results[0] = 0;
    for (Index i = 1; i < str.size(); ++i) {
        results[i] = results[i - 1] + (str[i - 1] == '1');
    }
    Bitvector<TypeParam> bv(str);
    Index numZeros = bv.rankZero(bv.sizeInBits() - 1);
    if (!bv.getBit(bv.sizeInBits() - 1)) {
        ++numZeros;
    }
    ASSERT_EQ(numZeros, bv.numZeros());
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.rankOne(i), results[i]) << i << " " << bv.sizeInBits();
        Index rank = bv.rankZero(i);
        if (rank < numZeros) {
            ASSERT_GE(bv.selectZero(rank), i) << bv.sizeInBits();
        } else {
            ASSERT_EQ(rank, numZeros);
        }
        if (i < numZeros) {
            ASSERT_EQ(bv.rankZero(bv.selectZero(i)), i) << bv.sizeInBits();
        }
    }
}

#ifdef ADS_HAS_CPP20
TYPED_TEST(ManyBitvecLayoutsTest, Constexpr) {
    using T = Bitvector<TypeParam>;
    static_assert(T("1").sizeInBits() == 1);
    static_assert(T("").sizeInBits() == 0);
    static_assert(!T("00011").getBit(2));
    static_assert(T("00011").getBit(3));
    static_assert(T("10101").rankZero(3) == 1);
    static_assert(T("10101").selectOne(1) == 2);
}
#endif
