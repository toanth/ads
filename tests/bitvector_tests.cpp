
#include "../include/bitvector/cache_efficient_rank_bitvec.hpp"
#include "../include/bitvector/efficient_rank_bitvec.hpp"
#include "gtest/gtest.h"

using namespace ads;

using SmallSuperBlocks = EfficientRankBitvec<1, (1 << 8) / 64, std::uint32_t>;
using NonsensicalSuperBlocks = EfficientRankBitvec<2, (1 << 14) / 64, std::uint64_t>;
using OneBlockPerSuperBlock = EfficientRankBitvec<64, 64>;
using LargeBlockCounts = EfficientRankBitvec<16, (1 << 16)>;

template<ADS_BITVEC_CONCEPT Bitvec>
class ManyBitvecsTest : public ::testing::Test {};

template<ADS_BITVEC_CONCEPT Bitvec>
class EfficientBitvecsTest : public ::testing::Test {};

template<typename = void>
using Bitvector = CacheEfficientRankBitvec; // TODO: Remove

using ManyBitvecs
        = ::testing::Types<CacheEfficientRankBitvec, EfficientRankBitvec<>, SmallSuperBlocks, NonsensicalSuperBlocks, OneBlockPerSuperBlock, LargeBlockCounts>;

using EfficientBitvecs = ::testing::Types<CacheEfficientRankBitvec, EfficientRankBitvec<>>;
TYPED_TEST_SUITE(ManyBitvecsTest, ManyBitvecs);
TYPED_TEST_SUITE(EfficientBitvecsTest, EfficientBitvecs);

TYPED_TEST(ManyBitvecsTest, ConstructionSizes) {
    {
        Allocation<> allocation(TypeParam::allocatedSizeInLimbsForLimbs(1));
        TypeParam bv = TypeParam::uninitializedForSize(64, allocation.memory());
        ASSERT_EQ(bv.numSuperblocks(), 1);
        ASSERT_EQ(bv.sizeInBits(), 64);
        ASSERT_EQ(bv.sizeInLimbs(), 1);
        ASSERT_EQ(bv.allocatedSizeInLimbs(), allocation.size());
        for (Index i = 0; i < 100'000; i += i / 10 + 1) {
            bv = TypeParam(i, 0);
            ASSERT_EQ(bv.sizeInBits(), i);
            ASSERT_EQ(bv.bitView().size(), bv.sizeInBits());
            ASSERT_GE(bv.sizeInLimbs(), (i + 63) / 64); // TODO: Make exact?
            ASSERT_EQ(bv.limbView().size(), bv.sizeInLimbs());
            ASSERT_EQ(bv.numSuperblocks(), (bv.sizeInLimbs() + bv.numLimbsInSuperblock() - 1) / bv.numLimbsInSuperblock());
        }
    }

    TypeParam bv("1");
    ASSERT_EQ(bv.numSuperblocks(), 1);
    ASSERT_EQ(bv.numBlocks(), 1);
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.sizeInLimbs(), 1);

    std::string text;
    text.reserve(10'000);
    for (Index i = 0; i < 10'000; ++i) {
        TypeParam bv(text);
        ASSERT_EQ(bv.sizeInBits(), i);
        ASSERT_GE(bv.sizeInLimbs(), (i + 63) / 64);
        ASSERT_EQ(bv.numSuperblocks(), (bv.sizeInLimbs() + bv.numLimbsInSuperblock() - 1) / bv.numLimbsInSuperblock());
        text.push_back('1');
    }
}

TYPED_TEST(ManyBitvecsTest, ConstructionElements) {
    TypeParam bv = TypeParam::uninitializedForSize(1);
    bv.setLimb(0, 1);
    bv.buildRankMetadata(0);
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.getLimb(0), 1);
    ASSERT_EQ(bv.getBit(0), 1);
    ASSERT_EQ(*bv.bitView().begin(), 1);
    ASSERT_EQ(bv.toString(), "1");
    ASSERT_EQ(bv, Bitvector<>("1"));
    ASSERT_GT(bv, Bitvector<>("0"));
    ASSERT_EQ(bv.numZeros(), 0);
    ASSERT_EQ(bv.numOnes(), 1);
    bv = TypeParam::uninitializedForSize(64);
    bv.setLimb(0, Elem(1) << 63);
    bv.buildRankMetadata(0);
    ASSERT_EQ(bv.sizeInBits(), 64);
    ASSERT_EQ(bv.getBit(0), 0);
    ASSERT_EQ(bv.getBit(63), 1);
    ASSERT_FALSE(*bv.bitView().begin());
    ASSERT_GT(bv, Bitvector<>("0"));
    ASSERT_EQ(bv.numZeros(), 63);
    ASSERT_EQ(bv.numOnes(), 1);
    std::vector<Elem> values{Elem(-1), Elem(-1)};
    bv = TypeParam(values, 70);
    ASSERT_EQ(bv.size(), 70);
    ASSERT_EQ(bv.numLimbs(), 2);
    ASSERT_EQ(bv.numOnes(), 70);
    ASSERT_EQ(bv.numZeros(), 0);
    bv = TypeParam(65, 0);
    ASSERT_EQ(bv.size(), 65);
    ASSERT_EQ(bv.numOnes(), 0);
    ASSERT_EQ(bv.numZeros(), 65);
}

TYPED_TEST(ManyBitvecsTest, ConstructionFromStringview) {
    TypeParam bv("01");
    ASSERT_EQ(bv.sizeInBits(), 2);
    ASSERT_EQ(bv.getBit(0), 0);
    ASSERT_EQ(bv.getBit(1), 1);
    ASSERT_EQ(bv.toString(), "01");
    ASSERT_EQ(bv.numZeros(), bv.numOnes());
    std::string s("101110011101");
    bv = TypeParam(s);
    ASSERT_EQ(bv.sizeInLimbs(), 1);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.getBit(i), s[i] == '1');
    }
    ASSERT_EQ(bv.numZeros(), 4);
    s = std::string(200, '1');
    s[63] = s[64 + 63] = s[64 * 2 + 63] = '0';
    Allocation<> allocation(TypeParam::allocatedSizeInLimbsForLimbs(roundUpDiv(s.size(), 64)));
    bv = TypeParam(s, 2, allocation.memory());
    ASSERT_EQ(bv.sizeInLimbs(), 4);
    ASSERT_EQ(bv.allocatedSizeInLimbs(), allocation.size());
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.getLimb(i), Elem(-1) / 2);
    }
    for (Index i = 0; i < 200 - 64 * 3; ++i) {
        ASSERT_TRUE(bv.getBit(64 * 3 + i));
    }
    bv = TypeParam(std::string(17, 'e'), 16);
    for (Index i = 0; i < 17 * 4; ++i) {
        ASSERT_EQ(bv.getBit(i), i % 4 != 3) << i;
    }
    bv = TypeParam(std::string(123, '2'), 4);
    for (Index i = 0; i < bv.size(); i += 2) {
        ASSERT_TRUE(bv.getBit(i)) << i;
        ASSERT_FALSE(bv.getBit(i + 1)) << i + 1;
    }
}

TYPED_TEST(ManyBitvecsTest, RankOnly1s) {
    TypeParam bv("111");
    ASSERT_EQ(bv.getLimb(0), Elem(0b111));
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
        ASSERT_EQ(bv.rankZero(i), 0) << i;
    }
    ASSERT_EQ(bv.numZeros(), 0);
}

TYPED_TEST(ManyBitvecsTest, RankSmall) {
    TypeParam bv("110011001100");
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

TYPED_TEST(ManyBitvecsTest, RankOneSuperblock) {
    std::string s(TypeParam::superblockSize(), '1');
    for (Index i = 0; i < s.size(); i += 2) {
        s[i] = '0';
    }
    TypeParam bv(s);
    ASSERT_EQ(bv.toString(), s);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 1) / 2) << i;
        ASSERT_EQ(bv.rankOne(i), i / 2) << i;
    }
    ASSERT_EQ(bv.numZeros(), bv.sizeInBits() / 2);
}

TYPED_TEST(ManyBitvecsTest, RankManySuperblocks) {
    std::string s(TypeParam::superblockSize() * 7 + 123456, '1');
    for (Index i = 0; i < s.size(); i += 3) {
        s[i] = '0';
    }
    TypeParam bv(s);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 2) / 3) << i;
        ASSERT_EQ(bv.rankOne(i), (2 * i) / 3) << i;
    }
}

TYPED_TEST(ManyBitvecsTest, SelectOnly1s) {
    TypeParam bv("1111");
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.selectOne(i), i) << i;
        //        ASSERT_EQ(bv.selectZero(i), -1) << i;
    }
    bv = TypeParam(std::string(200'000, '1'));
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
    }
}

TYPED_TEST(ManyBitvecsTest, SelectSmall) {
    std::string s("1100000000000000000000000000001");
    TypeParam bv(s);
    ASSERT_EQ(bv.selectOne(0), 0);
    ASSERT_EQ(bv.selectOne(1), 1);
    ASSERT_EQ(bv.selectOne(2), s.size() - 1);
    for (Index i = 0; i < s.size() - 3; ++i) {
        ASSERT_EQ(bv.selectZero(i), i + 2);
    }
    bv = TypeParam("0111");
    ASSERT_EQ(bv.selectOne(0), 1);
    ASSERT_EQ(bv.selectOne(1), 2);
    ASSERT_EQ(bv.selectZero(0), 0);
}

TYPED_TEST(ManyBitvecsTest, SelectLarge) {
    std::string s(12345, 'c');
    TypeParam bv(s, 16);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.selectZero(2 * i), 4 * i + 2);
        ASSERT_EQ(bv.selectZero(2 * i + 1), 4 * i + 3);
        ASSERT_EQ(bv.selectOne(2 * i), 4 * i);
        ASSERT_EQ(bv.selectOne(2 * i + 1), 4 * i + 1);
    }
}

TYPED_TEST(ManyBitvecsTest, Select2SuperBlocksPlus2) {
    TypeParam bv;
    bv = TypeParam(bv.superblockSize() * 2 + 2, 1);
    ASSERT_EQ(bv.rankOne(bv.sizeInBits() - 1), bv.sizeInLimbs());
    for (Index i = 0; i < bv.sizeInBits() - bv.sizeInLimbs(); ++i) {
        ASSERT_EQ(bv.selectZero(i), i + i / 63 + 1) << i << " " << bv.sizeInBits();
    }
    for (Index i = 0; i < bv.sizeInLimbs(); ++i) {
        ASSERT_EQ(bv.selectOne(i), i * 64) << i << " " << bv.sizeInBits();
    }
}


TYPED_TEST(ManyBitvecsTest, EmptyOrOneElem) {
    TypeParam bv = TypeParam::uninitializedForSize(0);
    ASSERT_EQ(bv.sizeInBits(), 0);
    ASSERT_EQ(bv.sizeInLimbs(), 0);
    bv = TypeParam(0, 123);
    ASSERT_EQ(bv.size(), 0);
    bv = TypeParam("1");
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.sizeInLimbs(), 1);
    ASSERT_EQ(bv.rankOne(0), 0);
    ASSERT_EQ(bv.rankZero(0), 0);
}

TYPED_TEST(EfficientBitvecsTest, PowerOfTwo) {
    for (Index i = 1; i <= (Elem(1) << 26); i *= 4) {
        TypeParam bv(i, Elem(0));
        for (Index j = 0; j < bv.numLimbs() - 1; ++j) {
            bv.setLimb(j, 0xaaaa'aaaa'aaaa'aaaaull);
        }
        for (Index j = (bv.numLimbs() - 1) * 64 + 1; j < bv.sizeInBits(); j += 2) {
            bv.setBit(j, true);
        }
        for (Index j = 0; j < bv.numSuperblocks(); ++j) {
            bv.buildRankMetadata(j);
        }
        for (Index j = 1; j < bv.sizeInBits(); j = j * 3 - 1) {
            ASSERT_EQ(bv.rankZero(j), (j + 1) / 2) << j << " " << i;
            ASSERT_EQ(bv.selectOne(j / 2), j / 2 * 2 + 1) << j << " " << i;
        }
        ASSERT_EQ(bv.numOnes(), bv.sizeInBits() / 2);
    }
}

TYPED_TEST(ManyBitvecsTest, Random) {
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
    TypeParam bv(str);
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

TYPED_TEST(ManyBitvecsTest, RandomLongRuns) {
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
    TypeParam bv(str);
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

// #ifdef ADS_HAS_CPP20
// TYPED_TEST(ManyBitvecsTest, Constexpr) {
//     static_assert(TypeParam("1").sizeInBits() == 1);
//     static_assert(TypeParam("").sizeInBits() == 0);
//     static_assert(!TypeParam("00011").getBit(2));
//     static_assert(TypeParam("00011").getBit(3));
//     static_assert(TypeParam("10101").rankZero(3) == 1);
//     static_assert(TypeParam("10101").selectOne(1) == 2);
//     static_assert(TypeParam(3, 0xf0f0'f0f0'f0f0'f0f0).selectOne(78) == 42);
// }
// #endif
