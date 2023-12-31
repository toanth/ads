
#include "../include/bitvector/cache_efficient_rank_bitvec.hpp"
#include "../include/bitvector/classical_bitvec.hpp"
#include "../include/bitvector/classical_rank_bitvec.hpp"
#include "../include/bitvector/recursive_bitvec.hpp"
// #include "../include/bitvector/select_bitvec.hpp"
#include "../include/bitvector/trivial_bitvec.hpp"
#include "gtest/gtest.h"

using namespace ads;

using SmallSuperBlocks = ClassicalRankBitvec<1, 8, std::uint32_t>;
using StrangeSuperBlocks = ClassicalRankBitvec<2, (1 << 14) / 64, std::uint64_t>;
using OneBlockPerSuperBlock = ClassicalRankBitvec<8, 8>;
using Classical = ClassicalBitvec<>;
using StrangeClassical = ClassicalBitvec<2, 64, 128 * 64, 256, U64, U32>;
using Trivial = TrivialBitvec<>;
using StrangeTrivial = TrivialBitvec<U64, ads::Operations::SELECT_ONLY>;
using VeryRecursive
        = RecursiveBitvec<RecursiveBitvec<RecursiveBitvec<RecursiveBitvec<CacheEfficientRankBitvec>>, SupportedSelects::ONE_ONLY, 64, U32>>;
using Strange = RecursiveBitvec<TrivialBitvec<U32, Operations::RANK_ONLY>, SupportedSelects::ONE_ONLY, 64 * 2>;
// using Select = SelectBitvec<>;

static_assert(IsNormalBitvec<VeryRecursive>);

template<ADS_BITVEC_CONCEPT Bitvec>
class AllBitvecsTest : public ::testing::Test {};

template<ADS_NORMAL_BITVEC_CONCEPT Bitvec>
class NormalBitvecsTest : public ::testing::Test {};

template<ADS_BITVEC_CONCEPT Bitvec>
class EfficientBitvecsTest : public ::testing::Test {};

template<ADS_BITVEC_CONCEPT Bitvec>
class EfficientNormalBitvecsTest : public ::testing::Test {};


template<typename = void>
using ReferenceBitvector = CacheEfficientRankBitvec;

/// \brief A collectio of bitvectors that try to achive high coverage
using AllBitvecs = ::testing::Types<CacheEfficientRankBitvec, ClassicalRankBitvec<>, OneBlockPerSuperBlock, Classical, StrangeClassical,
        Trivial, StrangeTrivial, EfficientBitvec<>, RecursiveBitvec<TrivialBitvec<>>, VeryRecursive, Strange /*, Select*/>;

/// \brief All bitvectors except TrivialBitvec, for which some operations aren't defined
using NormalBitvecs = ::testing::Types<CacheEfficientRankBitvec, ClassicalRankBitvec<>, SmallSuperBlocks,
        StrangeSuperBlocks, OneBlockPerSuperBlock, Classical, StrangeClassical, EfficientBitvec<>, Strange /*, Select*/>;

/// \brief Some of the faster bitvectors, which can be used for larger tests
using EfficientBitvecs = ::testing::Types<CacheEfficientRankBitvec, ClassicalRankBitvec<>, TrivialBitvec<U64>,
        EfficientBitvec<>, Classical, StrangeClassical, VeryRecursive, Strange /*, Select*/>;

using EfficientNormalBitvecs
        = ::testing::Types<CacheEfficientRankBitvec, ClassicalRankBitvec<>, Classical, StrangeClassical, EfficientBitvec<>, VeryRecursive /*, Select*/>;
TYPED_TEST_SUITE(AllBitvecsTest, AllBitvecs);
TYPED_TEST_SUITE(NormalBitvecsTest, NormalBitvecs);
TYPED_TEST_SUITE(EfficientBitvecsTest, EfficientBitvecs);
TYPED_TEST_SUITE(EfficientNormalBitvecsTest, EfficientNormalBitvecs);

TYPED_TEST(AllBitvecsTest, ConstructionSizes) {
    {
        Allocation<> allocation(TypeParam::allocatedSizeInBytesForBits(64));
        TypeParam bv = TypeParam::uninitializedForSize(64, allocation.memory());
        ASSERT_EQ(bv.size(), 64);
        if constexpr (IsSuperblockBitvec<TypeParam>) {
            ASSERT_EQ(bv.numSuperblocks(), 1);
            ASSERT_EQ(bv.template num<Group::CacheLine>(), 1);
        }
        ASSERT_EQ(bv.allocatedSizeInBytes(), allocation.sizeInBytes());
        ASSERT_EQ(bv.template allocatedSizeIn<Group::CacheLine>() * 8 * CACHELINE_SIZE_BYTES, bv.allocatedSizeInBits());
        for (Index i = 0; i < 100'000; i += i / 10 + 1) {
            bv = TypeParam(i, 0);
            ASSERT_EQ(bv.size(), i);
            ASSERT_EQ(bv.bitView().size(), bv.size());
            if constexpr (IsNormalBitvec<TypeParam>) {
                ASSERT_EQ(bv.numLimbs(), (i + 63) / 64);
                ASSERT_EQ(bv.limbView().size(), bv.numLimbs());
                if constexpr (IsSuperblockBitvec<TypeParam>) {
                    ASSERT_EQ(bv.numSuperblocks(), (bv.numLimbs() + bv.numLimbsInSuperblock() - 1) / bv.numLimbsInSuperblock());
                }
            }
        }
    }

    TypeParam bv("1");
    ASSERT_EQ(bv.size(), 1);
    if constexpr (IsNormalBitvec<TypeParam>) {
        ASSERT_EQ(bv.numLimbs(), 1);
        if constexpr (IsSuperblockBitvec<TypeParam>) {
            ASSERT_EQ(bv.numSuperblocks(), 1);
            ASSERT_EQ(bv.numBlocks(), 1);
        }
    }

    std::string text;
    text.reserve(10'000);
    for (Index i = 0; i < 10'000; ++i) {
        TypeParam bv(text);
        ASSERT_EQ(bv.size(), i);
        if constexpr (IsNormalBitvec<TypeParam>) {
            ASSERT_GE(bv.numLimbs(), (i + 63) / 64);
            if constexpr (IsSuperblockBitvec<TypeParam>) {
                ASSERT_EQ(bv.numSuperblocks(), (bv.numLimbs() + bv.numLimbsInSuperblock() - 1) / bv.numLimbsInSuperblock());
            }
        }
        text.push_back('1');
    }
}

TYPED_TEST(AllBitvecsTest, ConstructionElements) {
    TypeParam bv = TypeParam::uninitializedForSize(1);
    ASSERT_EQ(bv.size(), 1);
    if constexpr (IsNormalBitvec<TypeParam>) {
        bv.setLimb(0, 1);
        bv.buildMetadata();
        ASSERT_EQ(bv.getLimb(0), 1);
    } else {
        bv = TypeParam(1, Limb(1));
    }
    ASSERT_EQ(bv.getBit(0), 1);
    ASSERT_EQ(*bv.bitView().begin(), 1);
    ASSERT_EQ(bv.toString(), "1");
#ifdef ADS_HAS_CPP20
    ASSERT_EQ(bv, ReferenceBitvector<>("1"));
    ASSERT_GT(bv, ReferenceBitvector<>("0"));
#endif
    ASSERT_EQ(bv.numZeros(), 0);
    ASSERT_EQ(bv.numOnes(), 1);
    bv = TypeParam::uninitializedForSize(64);
    if constexpr (IsNormalBitvec<TypeParam>) {
        bv.setLimb(0, Limb(1) << 63);
        bv.buildMetadata();
    } else {
        ASSERT_EQ(bv.size(), 64);
        bv = TypeParam(64, Limb(1) << 63);
    }
    ASSERT_EQ(bv.size(), 64);
    ASSERT_EQ(bv.getBit(0), 0);
    ASSERT_EQ(bv.getBit(63), 1);
    ASSERT_FALSE(*bv.bitView().begin());
#ifdef ADS_HAS_CPP20
    ASSERT_GT(bv, ReferenceBitvector<>("0"));
#endif
    ASSERT_EQ(bv.numZeros(), 63);
    ASSERT_EQ(bv.numOnes(), 1);
    std::vector<Limb> values{Limb(-1), Limb(-1)};
    bv = TypeParam(values, 70);
    ASSERT_EQ(bv.size(), 70);
    if constexpr (IsNormalBitvec<TypeParam>) {
        ASSERT_EQ(bv.numLimbs(), 2);
    }
    ASSERT_EQ(bv.numOnes(), 70);
    ASSERT_EQ(bv.numZeros(), 0);
    bv = TypeParam(65, 0);
    ASSERT_EQ(bv.size(), 65);
    ASSERT_EQ(bv.numOnes(), 0);
    ASSERT_EQ(bv.numZeros(), 65);
    bv = TypeParam(64 * 23 + 1, Limb(9));
    ASSERT_EQ(bv.numOnes(), 23 * 2 + 1);
    ASSERT_TRUE(bv.getBit(bv.size() - 1));
}

TYPED_TEST(AllBitvecsTest, ConstructionFromStringview) {
    TypeParam bv("01");
    ASSERT_EQ(bv.size(), 2);
    ASSERT_EQ(bv.getBit(0), 0);
    ASSERT_EQ(bv.getBit(1), 1);
    ASSERT_EQ(bv.toString(), "01");
    ASSERT_EQ(bv.numZeros(), bv.numOnes());
    std::string s("101110011101");
    bv = TypeParam(s);
    ASSERT_EQ(bv.size(), ssize(s));
    for (Index i = 0; i < ssize(s); ++i) {
        ASSERT_EQ(bv.getBit(i), s[i] == '1');
    }
    ASSERT_EQ(bv.numZeros(), 4);
    s = std::string(200, '1');
    s[63] = s[64 + 63] = s[64 * 2 + 63] = '0';
    Allocation<> allocation(TypeParam::allocatedSizeInBytesForBits(ssize(s)));
    bv = TypeParam(s, 2, allocation.memory());
    ASSERT_EQ(bv.allocatedSizeInBytes(), allocation.sizeInBytes());
    if constexpr (IsNormalBitvec<TypeParam>) {
        ASSERT_EQ(bv.numLimbs(), 4);
        for (Index i = 0; i < 3; ++i) {
            ASSERT_EQ(bv.getLimb(i), Limb(-1) / 2);
        }
    }
    for (Index i = 0; i < 200 - 64 * 3; ++i) {
        ASSERT_TRUE(bv.getBit(64 * 3 + i)) << i;
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

TYPED_TEST(AllBitvecsTest, RankOnly1s) {
    TypeParam bv("111");
    if constexpr (IsNormalBitvec<TypeParam>) {
        ASSERT_EQ(bv.getLimb(0), Limb(0b111));
    }
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
        ASSERT_EQ(bv.rankZero(i), 0) << i;
    }
    ASSERT_EQ(bv.numZeros(), 0);
}

TYPED_TEST(AllBitvecsTest, RankSmall) {
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

TYPED_TEST(NormalBitvecsTest, RankOneSuperblock) {
    Index size = 1 << 16;
    if constexpr (IsSuperblockBitvec<TypeParam>) {
        size = TypeParam::superblockSize();
    }
    std::string s(size, '1');
    for (Index i = 0; i < ssize(s); i += 2) {
        s[i] = '0';
    }
    TypeParam bv(s);
    ASSERT_EQ(bv.toString(), s);
    for (Index i = 0; i < ssize(s); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 1) / 2) << i;
        ASSERT_EQ(bv.rankOne(i), i / 2) << i;
    }
    ASSERT_EQ(bv.numZeros(), bv.size() / 2);
}

TYPED_TEST(NormalBitvecsTest, RankManySuperblocks) {
    Index superblockSize = 1 << 16;
    if constexpr (IsSuperblockBitvec<TypeParam>) {
        superblockSize = TypeParam::superblockSize();
    }
    std::string s(superblockSize * 7 + 123456, '1');
    for (Index i = 0; i < ssize(s); i += 3) {
        s[i] = '0';
    }
    TypeParam bv(s);
    for (Index i = 0; i < ssize(s); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 2) / 3) << i;
        ASSERT_EQ(bv.rankOne(i), (2 * i) / 3) << i;
    }
}

TYPED_TEST(AllBitvecsTest, SelectSmall) {
    std::string s("1100000000000000000000000000001");
    TypeParam bv(s);
    ASSERT_EQ(bv.selectOne(0), 0);
    ASSERT_EQ(bv.selectOne(1), 1);
    ASSERT_EQ(bv.selectOne(2), ssize(s) - 1);
    for (Index i = 0; i < ssize(s) - 3; ++i) {
        ASSERT_EQ(bv.selectZero(i), i + 2);
    }
    bv = TypeParam("0111");
    ASSERT_EQ(bv.selectOne(0), 1);
    ASSERT_EQ(bv.selectOne(1), 2);
    ASSERT_EQ(bv.selectZero(0), 0);
}

TYPED_TEST(AllBitvecsTest, SelectLarge) {
    std::string s(12345, 'c');
    TypeParam bv(s, 16);
    for (Index i = 0; i < ssize(s); ++i) {
        ASSERT_EQ(bv.selectOne(2 * i), 4 * i) << i;
        ASSERT_EQ(bv.selectOne(2 * i + 1), 4 * i + 1) << i;
        ASSERT_EQ(bv.selectZero(2 * i), 4 * i + 2) << i;
        ASSERT_EQ(bv.selectZero(2 * i + 1), 4 * i + 3) << i;
    }
}

TYPED_TEST(NormalBitvecsTest, Select2SuperBlocksPlus2) {
    TypeParam bv;
    Index size = (1 << 17) + 42;
    if constexpr (IsSuperblockBitvec<TypeParam>) {
        size = bv.superblockSize() * 2 + 2;
    }
    bv = TypeParam(size, 1);
    ASSERT_EQ(bv.rankOne(bv.size() - 1), bv.numLimbs());
    for (Index i = 0; i < bv.size() - bv.numLimbs(); ++i) {
        ASSERT_EQ(bv.selectZero(i), i + i / 63 + 1) << i << " " << bv.size();
    }
    for (Index i = 0; i < bv.numLimbs(); ++i) {
        ASSERT_EQ(bv.selectOne(i), i * 64) << i << " " << bv.size();
    }
}

TYPED_TEST(AllBitvecsTest, SetBitsOneByOne) {
    Index size = 9876;
    TypeParam bv = TypeParam::uninitializedForSize(size);
    Index numOnes = 0;
    for (Index i = 0; i < size; ++i) {
        bool bit = ((i + 1) / 2) % 2;
        bv.setBit(i, bit);
        bv.buildMetadata();
        Index numZeros = i - numOnes;
        ASSERT_EQ(bv.size(), size);
        ASSERT_EQ(bv.getBit(i), bit) << i;
        ASSERT_EQ(bv.rankOne(i), numOnes) << i;
        ASSERT_EQ(bv.rankZero(i), i - numOnes) << i;
        if (numOnes + bit > 0) {
            if (bit) {
                ASSERT_EQ(bv.selectOne(numOnes), i) << i << ", " << numOnes;
            } else {
                ASSERT_LT(bv.selectOne(numOnes - 1), i) << i << ", " << numOnes;
            }
        }
        if (numZeros + 1 - bit > 0) {
            if (!bit) {
                ASSERT_EQ(bv.selectZero(numZeros), i);
            } else {
                ASSERT_LT(bv.selectZero(numZeros - 1), i);
            }
        }
        numOnes += bit;
    }
    ASSERT_EQ(bv.numOnes(), numOnes);
    ASSERT_EQ(bv.numZeros(), size - numOnes);
}

TYPED_TEST(AllBitvecsTest, Only1s) {
    TypeParam bv("1111");
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
        ASSERT_EQ(bv.selectOne(i), i) << i;
    }
    bv = TypeParam(std::string(200'000, '1'));
    ASSERT_EQ(bv.numZeros(), 0);
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
        ASSERT_EQ(bv.selectOne(i), i) << i;
    }
    bv = TypeParam::uninitializedForSize(65538);
    for (Index i = 0; i < 65538; ++i) {
        bv.setBit(i);
    }
    bv.buildMetadata();
    ASSERT_EQ(bv.numOnes(), bv.size());
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankOne(i), i);
        ASSERT_EQ(bv.selectOne(i), i);
    }
}

TYPED_TEST(AllBitvecsTest, Only0s) {
    TypeParam bv("00000");
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), i) << i;
        ASSERT_EQ(bv.selectZero(i), i) << i;
    }
    bv = TypeParam(std::string(234'567, '0'));
    ASSERT_EQ(bv.numOnes(), 0);
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), i) << i;
        ASSERT_EQ(bv.selectZero(i), i) << i;
    }
    bv = TypeParam::uninitializedForSize(65537);
    for (Index i = 0; i < 65537; ++i) {
        bv.setBit(i, false);
    }
    bv.buildMetadata();
    ASSERT_EQ(bv.numZeros(), bv.size());
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), i);
        ASSERT_EQ(bv.selectZero(i), i);
    }
}

TYPED_TEST(AllBitvecsTest, AlternatingOneZero) {
    TypeParam bv(90'001, 0xaaaa'aaaa'aaaa'aaaaull);
    ASSERT_EQ(bv.numZeros(), bv.numOnes() + 1);
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankOne(i), i / 2) << i;
        ASSERT_EQ(bv.rankZero(i), (i + 1) / 2) << i;
    }
    for (Index i = 0; i < bv.numOnes(); ++i) {
        ASSERT_EQ(bv.selectOne(i), 2 * i + 1) << i;
        ASSERT_EQ(bv.selectZero(i), 2 * i) << i;
    }
    ASSERT_EQ(bv.selectZero(bv.size() / 2), bv.size() - 1);
}

TYPED_TEST(AllBitvecsTest, IncreasingRunLengths) {
    TypeParam bv = TypeParam::uninitializedForSize(520 * 519 / 2);
    for (Index i = 0, current = 0; i < 520; ++i) {
        for (Index j = 0; j < i; ++j, ++current) {
            bv.setBit(current, i % 2);
        }
    }
    bv.buildMetadata();

    Index current = 0;
    Index numOnes = 0;
    for (Index i = 0; i < 520; ++i) {
        for (Index j = 0; j < i; ++j, ++current, numOnes += i % 2) {
            ASSERT_EQ(bv.getBit(current), i % 2) << i << " " << j << " " << current << " " << numOnes;
            ASSERT_EQ(bv.rankOne(current), numOnes) << i << " " << j << " " << current << " " << numOnes;
            ASSERT_EQ(bv.rankZero(current), current - numOnes) << i << " " << j << " " << current << " " << numOnes;
            if (i % 2) {
                ASSERT_EQ(bv.selectOne(numOnes), current) << i << " " << j << " " << current << " " << numOnes;
            } else {
                ASSERT_EQ(bv.selectZero(current - numOnes), current) << i << " " << j << " " << current << " " << numOnes;
            }
        }
    }
    ASSERT_EQ(current, bv.size());
    ASSERT_EQ(bv.numOnes(), numOnes);
    ASSERT_EQ(bv.numZeros(), bv.size() - numOnes);
}

TYPED_TEST(NormalBitvecsTest, ExponentialOneIndices) {
    TypeParam bv = TypeParam(98765, 0);
    Index numOnes = intLog2(bv.size());
    for (Index i = 0; i < numOnes; ++i) {
        bv.setBit(Elem(1) << i);
    }
    bv.buildMetadata();
    for (Index i = 0; i < numOnes; ++i) {
        ASSERT_EQ(bv.selectOne(i), Elem(1) << i);
        ASSERT_LE(bv.selectZero(i), i + 2 + intLog2(i + 1)) << i;
        ASSERT_EQ(bv.rankOne((Elem(1) << i) + 1), i + 1);
        if (i > 1) {
            ASSERT_EQ(bv.rankOne((Elem(1) << i) - 1), i);
            ASSERT_EQ(bv.rankOne((Elem(1) << i) + (Elem(1) << (i - 1))), i + 1);
        }
    }
}

TYPED_TEST(AllBitvecsTest, EmptyOrOneElem) {
    TypeParam bv = TypeParam::uninitializedForSize(0);
    ASSERT_EQ(bv.size(), 0);
    if constexpr (IsNormalBitvec<TypeParam>) {
        ASSERT_EQ(bv.numLimbs(), 0);
    }
    bv = TypeParam(0, 123);
    ASSERT_EQ(bv.size(), 0);
    bv = TypeParam("1");
    ASSERT_EQ(bv.size(), 1);
    if constexpr (IsNormalBitvec<TypeParam>) {
        ASSERT_EQ(bv.numLimbs(), 1);
    }
    ASSERT_EQ(bv.rankOne(0), 0);
    ASSERT_EQ(bv.rankZero(0), 0);
}

TYPED_TEST(EfficientBitvecsTest, PowerOfTwo) {
    for (Index i = 1; i <= (Index(1) << 26); i *= 4) {
        TypeParam bv(i, Limb(0));
        if constexpr (IsSuperblockBitvec<TypeParam>) {
            for (Index j = 0; j < bv.numLimbs() - 1; ++j) {
                bv.setLimb(j, 0xaaaa'aaaa'aaaa'aaaaull);
            }
            for (Index j = (bv.numLimbs() - 1) * 64 + 1; j < bv.size(); j += 2) {
                bv.setBit(j, true);
            }
            for (Index j = 0; j < bv.numSuperblocks(); ++j) {
                bv.buildRankMetadata(j);
            }
        } else {
            bv = TypeParam(i, Limb(0xaaaa'aaaa'aaaa'aaaa));
        }
        for (Index j = 1; j < bv.size(); j = j * 3 - 1) {
            ASSERT_EQ(bv.rankZero(j), (j + 1) / 2) << j << " " << i;
            ASSERT_EQ(bv.selectOne(j / 2), j / 2 * 2 + 1) << j << " " << i;
        }
        ASSERT_EQ(bv.numOnes(), bv.size() / 2) << i;
    }
}

TYPED_TEST(AllBitvecsTest, RandomSmallish) {
    auto engine = createRandomEngine();
    auto dist = std::uniform_int_distribution<Index>(65530, 65550);
    std::string str(dist(engine), '0');
    for (char& c : str) {
        if (dist(engine) & 1) {
            c = '1';
        }
    }
    std::vector<std::uint32_t> results(ssize(str));
    results[0] = 0;
    for (Index i = 1; i < ssize(str); ++i) {
        results[i] = results[i - 1] + (str[i - 1] == '0');
    }
    TypeParam bv(str);
    Index numOnes = bv.rankOne(bv.size() - 1);
    if (bv.getBit(bv.size() - 1)) {
        ++numOnes;
    }
    ASSERT_EQ(numOnes, bv.numOnes());
    ASSERT_EQ(bv.numZeros(), results.back() + (bv.getBit(bv.size() - 1) ? 0 : 1));
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), results[i]) << i << " " << bv.size();
        Index rank = bv.rankOne(i);
        ASSERT_EQ(rank, bv.template inefficientRank<true>(i)) << i;
        if (rank < numOnes) {
            Index j = bv.selectOne(rank);
            ASSERT_GE(j, i) << bv.size() << i;
            ASSERT_EQ(j, bv.template inefficientSelect<true>(rank)) << i;
        } else {
            ASSERT_EQ(rank, numOnes) << i;
        }
        if (i < numOnes) {
            ASSERT_EQ(bv.rankOne(bv.selectOne(i)), i) << bv.size();
        }
    }
}

TYPED_TEST(AllBitvecsTest, Random) {
    auto engine = createRandomEngine();
    auto dist = std::uniform_int_distribution<Index>(0, 400'000);
    std::string str(dist(engine), '0');
    for (char& c : str) {
        if (dist(engine) & 1) {
            c = '1';
        }
    }
    std::vector<std::uint32_t> results(ssize(str));
    results[0] = 0;
    for (Index i = 1; i < ssize(str); ++i) {
        results[i] = results[i - 1] + (str[i - 1] == '0');
    }
    TypeParam bv(str);
    Index numOnes = bv.rankOne(bv.size() - 1);
    if (bv.getBit(bv.size() - 1)) {
        ++numOnes;
    }
    ASSERT_EQ(numOnes, bv.numOnes());
    ASSERT_EQ(bv.numZeros(), results.back() + (bv.getBit(bv.size() - 1) ? 0 : 1));
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), results[i]) << i << " " << bv.size();
        Index rank = bv.rankOne(i);
        if (rank < numOnes) {
            ASSERT_GE(bv.selectOne(rank), i) << bv.size();
        } else {
            ASSERT_EQ(rank, numOnes);
        }
        if (i < numOnes) {
            ASSERT_EQ(bv.rankOne(bv.selectOne(i)), i) << bv.size();
        }
    }
}

TYPED_TEST(AllBitvecsTest, RandomLongRuns) {
    auto engine = createRandomEngine();
    auto dist = std::uniform_real_distribution<double>(-16.0, 16.0);
    Index numOnes = 0;
    std::string str;
    while (ssize(str) < 300'000) {
        double randomVal = dist(engine);
        Index len = Index(std::max(0.0, std::pow(2.0, std::abs(randomVal))));
        char c = randomVal > 0.0 ? '1' : '0';
        str.append(len, c);
        numOnes += (randomVal > 0.0 ? len : 0);
    }
    std::vector<std::uint32_t> rankOnes(ssize(str) + 1);
    rankOnes[0] = 0;
    for (Index i = 1; i < ssize(str) + 1; ++i) {
        rankOnes[i] = rankOnes[i - 1] + (str[i - 1] == '1');
    }
    ASSERT_EQ(rankOnes[ssize(str)], numOnes);
    TypeParam bv(str);
    Index numZeros = bv.rankZero(bv.size() - 1);
    if (!bv.getBit(bv.size() - 1)) {
        ++numZeros;
    }
    ASSERT_EQ(numZeros, bv.numZeros());
    ASSERT_EQ(bv.numOnes(), numOnes);
    ASSERT_EQ(bv.numOnes() + bv.numZeros(), bv.size());
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.rankOne(i), rankOnes[i]) << i << " " << bv.size();
        Index rank = bv.rankZero(i);
        if (rank < numZeros) {
            ASSERT_GE(bv.selectZero(rank), i) << bv.size();
        } else {
            ASSERT_EQ(rank, numZeros);
        }
        if (i < numZeros) {
            ASSERT_EQ(bv.rankZero(bv.selectZero(i)), i) << bv.size();
        }
    }
}

TYPED_TEST(EfficientNormalBitvecsTest, RandomFewOneClusters) {
    auto engine = createRandomEngine();
    auto sizes = std::uniform_int_distribution<Index>(1 << 28, Elem(1) << 31);
    auto limbs = std::uniform_int_distribution<Limb>();
    TypeParam bv = TypeParam::uninitializedForSize(sizes(engine));
    Index prob = sizes(engine) % 1233 + 100;
    for (Index i = 0; i < bv.numLimbs(); ++i) {
        Limb limb = limbs(engine);
        if (limb % prob != 0) [[likely]] {
            limb = 0;
        }
        bv.setLimb(i, limb);
    }
    bv.buildMetadata();
    ASSERT_GE(bv.numZeros(), bv.numOnes()); // vanishingly small chance of failure
    for (Index i = 0; i < bv.numOnes(); ++i) {
        ASSERT_EQ(bv.rankOne(bv.selectOne(i)), i);
        if (i % prob < 5) {
            ASSERT_EQ(bv.template inefficientSelect<true>(i), bv.selectOne(i));
        }
        ASSERT_LE(bv.selectZero(bv.rankOne(i)), i / 10 + 100); // unlikely to fail
    }
}

TYPED_TEST(AllBitvecsTest, InefficientOps) {
    TypeParam bv("1111000111");
    ASSERT_EQ(bv.template inefficientRank<true>(2), 2);
    ASSERT_EQ(bv.template inefficientRank<true>(6), 4);
    ASSERT_EQ(bv.template inefficientRank<false>(4), 0);
    ASSERT_EQ(bv.template inefficientSelect<true>(1), 1);
    ASSERT_EQ(bv.template inefficientSelect<false>(1), 5);
    bv = TypeParam(65540, 0x5555'5555'5555'5555ull);
    auto select0 = bv.template selectView<false>();
    auto select1 = bv.template selectView<true>();
    auto rank0 = bv.template rankView<false>();
    auto rank1 = bv.template rankView<true>();
    ASSERT_EQ(rank0.size(), bv.size());
    ASSERT_EQ(rank1.size(), bv.size());
    ASSERT_EQ(select0.size() * 2, bv.size());
    ASSERT_EQ(select1.size() * 2, bv.size());
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(rank0[i], i / 2) << i;
        ASSERT_EQ(rank1[i], (i + 1) / 2) << i;
    }
    ASSERT_EQ(bv.selectZero(0), 1);
    ASSERT_EQ(selectFunc<false>(bv, 0), 1);
    ASSERT_EQ(select0[0], 1);
    for (Index i = 0; i < bv.size() / 2; ++i) {
        ASSERT_EQ(select0[i], i * 2 + 1) << i;
        ASSERT_EQ(select1[i], i * 2) << i;
    }
    for (Index i = 0; i < bv.size(); ++i) {
        ASSERT_EQ(bv.getBit(i), (i + 1) % 2);
        ASSERT_EQ(bv.template inefficientRank<false>(i), i / 2) << i;
        ASSERT_EQ(bv.template inefficientRank<true>(i), (i + 1) / 2) << i;
    }
    for (Index i = 0; i * 2 < bv.size(); ++i) {
        ASSERT_EQ(bv.template inefficientSelect<false>(i), i * 2 + 1) << i;
        ASSERT_EQ(bv.template inefficientSelect<true>(i), i * 2) << i;
    }
}


template<ADS_BITVEC_CONCEPT Bitvec>
static ADS_CONSTEVAL Index constexprBitvecTests() noexcept {
    Bitvec bv("10010");
    assert(bv.rankZero(3) == 2);
    assert(bv.selectOne(1) == 3);
    bv = Bitvec::uninitializedForSize(123);
    for (Index i = 0; i < bv.size(); ++i) {
        bv.setBit(i, i % 3 == 0);
    }
    bv.buildMetadata();
    for (Index i = 0; i < bv.size(); ++i) {
        assert(bv.rankOne(i) == (i + 2) / 3);
    }
    for (Index i = 0; i < bv.numOnes(); ++i) {
        assert(bv.selectOne(i) == i * 3);
        assert(bv.selectZero(2 * i) == (3 * i + 1));
        assert(bv.selectZero(2 * i + 1) == (3 * i + 2));
    }
    if constexpr (IsSuperblockBitvec<Bitvec>) {
        bv = Bitvec(Bitvec::superblockSize() * 2 + 1, Limb(5));
        assert(bv.numOnes() == bv.numLimbs() * 2 - 1);
        assert(bv.selectOne(bv.numLimbsInSuperblock() * 2) == bv.superblockSize());
    } else if constexpr (IsNormalBitvec<Bitvec>) {
        bv = Bitvec(Bitvec::blockSize() * 2 + 1, Limb(5));
        assert(bv.numOnes() == bv.numLimbs() * 2 - 1);
        assert(bv.selectOne(bv.numLimbsInBlock() * 2) == bv.blockSize());
    } else {
        bv = Bitvec(64 * 17 + 1, Limb(5));
        assert(bv.selectOne(21) == 10 * 64 + 2);
    }
    assert(bv.numOnes() * 62 > bv.numZeros() * 2);
    assert((bv.numOnes() - 1) * 62 == bv.numZeros() * 2);
    assert(bv.getBit(bv.size() - 1));
    assert(bv.rankZero(bv.size() - 1) == bv.numZeros());
    return 0;
}

#ifdef ADS_HAS_CPP20
TYPED_TEST(AllBitvecsTest, Constexpr) {
    static_assert(TypeParam("1").size() == 1);
    static_assert(TypeParam("").size() == 0);
    static_assert(!TypeParam("00011").getBit(2));
    static_assert(TypeParam("00011").getBit(3));
    static_assert(TypeParam("10101").rankZero(3) == 1);
    static_assert(TypeParam("10101").selectOne(1) == 2);
    static_assert(TypeParam(3 * 64 - 7, 0xf0f0'f0f0'f0f0'f0f0).selectOne(42) == 86);
    static_assert(constexprBitvecTests<TypeParam>() == 0);
}
#endif
