#include "../include/bitvector.hpp"
#include "gtest/gtest.h"

using namespace ads;

// using TestLayout = CacheEfficientLayout;
using TestLayout = SimpleLayout<>;

TEST(BitvectorConstruction, Sizes) {
    {
        Bitvector<TestLayout> bv(64);
        ASSERT_EQ(bv.numSuperBlocks(), 1);
        ASSERT_EQ(bv.sizeInBits(), 64);
        ASSERT_EQ(bv.sizeInElems(), 1);
        for (Index i = 0; i < 100'000; ++i) {
            bv = Bitvector<TestLayout>(i);
            ASSERT_EQ(bv.sizeInBits(), i);
            ASSERT_EQ(bv.bitView().size(), bv.sizeInBits());
            ASSERT_GE(bv.sizeInElems(), (i + 63) / 64); // TODO: Make exact?
            ASSERT_EQ(bv.elemView().size(), bv.sizeInElems());
            ASSERT_EQ(bv.numSuperBlocks(), (bv.sizeInElems() + bv.numElemsInSuperBlock() - 1) / bv.numElemsInSuperBlock());
        }
    }

    Bitvector<TestLayout> bv("1");
    ASSERT_EQ(bv.numSuperBlocks(), 1);
    ASSERT_EQ(bv.numBlocks(), 1);
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.sizeInElems(), 1);

    std::string text;
    text.reserve(10'000);
    for (Index i = 0; i < 10'000; ++i) {
        Bitvector<TestLayout> bv(text);
        ASSERT_EQ(bv.sizeInBits(), i);
        ASSERT_GE(bv.sizeInElems(), (i + 63) / 64);
        ASSERT_EQ(bv.numSuperBlocks(), (bv.sizeInElems() + bv.numElemsInSuperBlock() - 1) / bv.numElemsInSuperBlock());
        text.push_back('1');
    }
}

TEST(BitvectorConstruction, Elements) {
    Bitvector<TestLayout> bv(1);
    bv.setElem(0, 1);
    bv.buildRankMetadata(0);
    ASSERT_EQ(bv.sizeInBits(), 1);
    ASSERT_EQ(bv.getElem(0), 1);
    ASSERT_EQ(bv.getBit(0), 1);
    ASSERT_EQ(*bv.bitView().begin(), 1);
    ASSERT_EQ(bv.toString(), "1");
    ASSERT_EQ(bv, Bitvector<>("1"));
    ASSERT_GT(bv, Bitvector<>("0"));
    bv = Bitvector(64);
    bv.setElem(0, Elem(1) << 63);
    bv.buildRankMetadata(0);
    ASSERT_EQ(bv.sizeInBits(), 64);
    ASSERT_EQ(bv.getBit(0), 0);
    ASSERT_EQ(bv.getBit(63), 1);
    ASSERT_FALSE(*bv.bitView().begin());
    ASSERT_GT(bv, Bitvector<>("0"));
}

TEST(BitvectorConstruction, FromStringview) {
    Bitvector<TestLayout> bv("01");
    ASSERT_EQ(bv.sizeInBits(), 2);
    ASSERT_EQ(bv.getBit(0), 0);
    ASSERT_EQ(bv.getBit(1), 1);
    for (Index i = 2; i < 64; ++i) {
        ASSERT_EQ(bv.getBit(i), 0);
    }
    ASSERT_EQ(bv.toString(), "01");
    std::string s("101110011101");
    bv = Bitvector<TestLayout>(s);
    ASSERT_EQ(bv.sizeInElems(), 1);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.getBit(i), s[i] == '1');
    }
    s = std::string(200, '1');
    s[63] = s[64 + 63] = s[64 * 2 + 63] = '0';
    bv = Bitvector<TestLayout>(s);
    ASSERT_EQ(bv.sizeInElems(), 4);
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(bv.getElem(i), Elem(-1) / 2);
    }
    for (Index i = 0; i < 64; ++i) {
        ASSERT_EQ(bv.getBit(64 * 3 + i), i < 200 - 64 * 3);
    }
    bv = Bitvector<TestLayout>(std::string(11, 'e'), 16);
    for (Index i = 0; i < 11 * 4; ++i) {
        ASSERT_EQ(bv.getBit(i), i % 4 != 3) << i;
    }
}

TEST(BitvecRank, Only1s) {
    Bitvector<TestLayout> bv("111");
    ASSERT_EQ(bv.getElem(0), Elem(0b111));
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
    std::string s(bv.superBlockSize(), '1');
    for (Index i = 0; i < s.size(); i += 2) {
        s[i] = '0';
    }
    bv = Bitvector<TestLayout>(s);
    ASSERT_EQ(bv.toString(), s);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 1) / 2);
        ASSERT_EQ(bv.rankOne(i), i / 2);
    }
}

TEST(BitvecRank, ManySuperblocks) {
    Bitvector<TestLayout> bv(0);
    std::string s(bv.superBlockSize() * 7 + 123456, '1');
    for (Index i = 0; i < s.size(); i += 3) {
        s[i] = '0';
    }
    bv = Bitvector<TestLayout>(s);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.rankZero(i), (i + 2) / 3) << i;
        ASSERT_EQ(bv.rankOne(i), (2 * i) / 3) << i;
    }
}

TEST(BitvecSelect, Only1s) {
    Bitvector<TestLayout> bv("1111");
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.selectOne(i), i) << i;
        //        ASSERT_EQ(bv.selectZero(i), -1) << i;
    }
    bv = Bitvector<TestLayout>(std::string(200'000, '1'));
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.rankOne(i), i) << i;
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
    bv = Bitvector<TestLayout>("0111");
    ASSERT_EQ(bv.selectOne(0), 1);
    ASSERT_EQ(bv.selectOne(1), 2);
    ASSERT_EQ(bv.selectZero(0), 0);
}

TEST(BitvecSelect, Large) {
    std::string s(12345, 'c');
    Bitvector<TestLayout> bv(s, 16);
    for (Index i = 0; i < s.size(); ++i) {
        ASSERT_EQ(bv.selectZero(2 * i), 4 * i + 2);
        ASSERT_EQ(bv.selectZero(2 * i + 1), 4 * i + 3);
        ASSERT_EQ(bv.selectOne(2 * i), 4 * i);
        ASSERT_EQ(bv.selectOne(2 * i + 1), 4 * i + 1);
    }
}

TEST(BitvecSelect, 2SuperBlocksPlus2) {
    Bitvector<TestLayout> bv;
    bv = Bitvector<TestLayout>(bv.superBlockSize() * 2 + 2, 1);
    ASSERT_EQ(bv.rankOne(bv.sizeInBits() - 1), bv.sizeInElems());
    for (Index i = 0; i < bv.sizeInBits() - bv.sizeInElems(); ++i) {
        ASSERT_EQ(bv.selectZero(i), i + i / 63 + 1) << i << " " << bv.sizeInBits();
    }
    for (Index i = 0; i < bv.sizeInElems(); ++i) {
        ASSERT_EQ(bv.selectOne(i), i * 64) << i << " " << bv.sizeInBits();
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

TEST(Bitvector, PowerOfTwo) {
    for (Index i = 1; i <= (Elem(1) << 26); i *= 4) {
        Bitvector<TestLayout> bv(i, 0);
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
    }
}

TEST(Bitvector, Random) {
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
    Bitvector<TestLayout> bv(str);
    Index numOnes = bv.rankOne(bv.sizeInBits() - 1);
    if (bv.getBit(bv.sizeInBits() - 1)) {
        ++numOnes;
    }
    for (Index i = 0; i < bv.sizeInBits(); ++i) {
        ASSERT_EQ(bv.rankZero(i), results[i]) << i << " " << bv.sizeInBits();
        Index rank = bv.rankOne(i);
        if (rank < numOnes) {
            ASSERT_GE(bv.selectOne(rank), i) << bv.sizeInBits();
        } else {
            ASSERT_EQ(rank, numOnes);
            //            ASSERT_EQ(bv.selectOne(rank), -1);
        }
        if (i < numOnes) {
            ASSERT_EQ(bv.rankOne(bv.selectOne(i)), i) << bv.sizeInBits();
        }
    }
}
