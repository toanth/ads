

#include "../include/succinct_rmq.hpp"
#include "gtest/gtest.h"

using namespace ads;

using Bitvec = EfficientSelectBitvec<>;

SuccinctRMQ<> randomRmq(Index length) {
    auto engine = createRandomEngine();
    auto dist = std::uniform_int_distribution<Index>();
    std::vector<Index> vec(length);
    for (Index& v : vec) {
        v = dist(engine);
    }
    return SuccinctRMQ<>(vec);
}

void testRmmTree(const RangeMinMaxTree<>& tree) {
    ASSERT_EQ(tree.size - 1, tree.leafIdxInArr(tree.numLeaves - 1));
    Index treeDepth = intLog2(tree.size - 1);
    for (Index i = 1; i < tree.leafIdxInArr(0); ++i) {
        Index leftChild = 2 * i;
        Index rightChild = 2 * i + 1;
        Index childDepth = intLog2(Elem(rightChild));
        bool rightChildExists = ((tree.size - 1) >> (treeDepth - childDepth)) >= rightChild;
        bool leftChildExists = ((tree.size - 1) >> (treeDepth - childDepth)) >= leftChild;
        if (rightChildExists) {
            ASSERT_EQ(tree.rmmArr[i], std::min(tree.rmmArr[i * 2], tree.rmmArr[i * 2 + 1]));
        } else if (leftChildExists) {
            assert(2 * i <= tree.leafIdxInArr(tree.numLeaves - 1));
            ASSERT_EQ(tree.rmmArr[i], tree.rmmArr[2 * i]);
        }
    }
}

TEST(RmmTree, Small) {
    Allocation allocation(Bitvec::allocatedSizeInLimbsForLimbs(1));
    RangeMinMaxTree<> tree(Bitvec("11110000"), allocation.memory());
    ASSERT_EQ(tree.size, 2);
    ASSERT_EQ(tree.findMinInBlock(0, 8).pos, 7);
    ASSERT_EQ(tree.findMinInBlock(0, 8).minExcess, 0);
    ASSERT_EQ(tree.findMinInBlock(0, 4).pos, 0);
    ASSERT_EQ(tree.findMinInBlock(4, 7).pos, 6);
    ASSERT_EQ(tree.bitvecRmq(0, 1), 0);
    ASSERT_EQ(tree.bitvecRmq(0, 2), 0);
    testRmmTree(tree);
}

TEST(RmmTree, 1BlockDescending) {
    Index blockSize = RangeMinMaxTree<>::blockSize;
    std::string str(blockSize, '1');
    for (Index i = str.size() / 2; i < str.size(); ++i) {
        str[i] = '0';
    }
    Allocation allocation(Bitvec::allocatedSizeInLimbsForBits(str.size()));
    RangeMinMaxTree<> tree((Bitvec(str)), allocation.memory());
    ASSERT_EQ(tree.size, 2);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(0)], 0);
    ASSERT_EQ(tree.findMinInTree(0, 1).minExcess, tree.findMinInBlock(0, blockSize).minExcess);
    for (Index i = 0; i < str.size(); ++i) {
        for (Index j = i + 1; j < str.size(); ++j) {
            Index e = tree.findMinInBlock(i, j).minExcess;
            if (j <= blockSize / 2) {
                ASSERT_EQ(e, i + 1) << i << " " << j;
            } else if (i >= blockSize / 2) {
                ASSERT_EQ(e, blockSize / 2 - j + blockSize / 2) << i << " " << j;
            } else {
                ASSERT_EQ(e, std::min(Index(i + 1), blockSize / 2 - i - j + blockSize / 2 + i)) << i << " " << j;
            }
        }
    }
    for (Index j = 0; j < str.size() / 2 - 1; ++j) {
        for (Index i = 0; i < j; ++i) {
            ASSERT_EQ(tree.bitvecRmq(i, j), i);
        }
    }
    testRmmTree(tree);
}

TEST(RmmTree, 1Block2BitsRandom) {
    SuccinctRMQ _ = randomRmq(RangeMinMaxTree<>::blockSize / 2);
    const RangeMinMaxTree<>& tree = _.getTree();
    ASSERT_EQ(tree.numLeaves, 2);
    ASSERT_EQ(tree.getBitvector().getBit(0), openParen);
    ASSERT_EQ(tree.getBitvector().getBit(tree.getBitvector().sizeInBits() - 1), closeParen);
    auto m = tree.findMinInBlock(0, RangeMinMaxTree<>::blockSize);
    ASSERT_EQ(m.pos, tree.findMinInTree(0, 1).pos);
    ASSERT_EQ(m.minExcess, tree.findMinInTree(0, 1).minExcess);
    if (tree.getBitvector().getBit(tree.getBitvector().sizeInBits() - 2) == openParen) {
        ASSERT_EQ(m.minExcess, 1);
    }
    ASSERT_LE(m.minExcess, 2);
    ASSERT_GE(m.minExcess, 1);
    ASSERT_EQ(tree.bitvecRmq(0, tree.getBitvector().sizeInBits()), tree.getBitvector().sizeInBits() - 1);
    testRmmTree(tree);
}

TEST(RmmTree, 2Blocks) {
    Index blockSize = RangeMinMaxTree<>::blockSize;
    std::string str(blockSize + 10, '1');
    for (Index i = str.size() / 2; i < str.size(); ++i) {
        str[i] = '0';
    }
    Allocation allocation(Bitvec::allocatedSizeInLimbsForBits(str.size()));
    RangeMinMaxTree<> tree((Bitvec(str)), allocation.memory());
    ASSERT_EQ(tree.bitvecRmq(blockSize - 4, str.size()), str.size() - 1);
    testRmmTree(tree);
    str = std::string(blockSize * 2, '1');
    for (Index i = 2; i < str.size(); i += 2) {
        str[i] = '0';
    }
    str.back() = '0';
    allocation = Allocation(Bitvec::allocatedSizeInLimbsForBits(str.size()));
    tree = RangeMinMaxTree<>((Bitvec(str)), allocation.memory());
    ASSERT_EQ(tree.bitvecRmq(0, blockSize), 0);
    ASSERT_EQ(tree.bitvecRmq(blockSize - 1, 2 * blockSize - 1), blockSize);
    ASSERT_EQ(tree.bitvecRmq(blockSize, 2 * blockSize), 2 * blockSize - 1);
    ASSERT_EQ(tree.bitvecRmq(2, 4), 2);
    ASSERT_EQ(tree.bitvecRmq(blockSize - 2, blockSize + 2), blockSize - 2);
    testRmmTree(tree);
}


TEST(RmmTree, 10BlocksAscending) {
    Index blockSize = RangeMinMaxTree<>::blockSize;
    std::string str(blockSize * 10, '1');
    for (Index i = str.size() / 2; i < str.size(); ++i) {
        str[i] = '0';
    }
    Allocation allocation(Bitvec::allocatedSizeInLimbsForBits(str.size()));
    RangeMinMaxTree<> tree((Bitvec(str)), allocation.memory());
    ASSERT_EQ(tree.bv.getBit(5 * blockSize - 1), openParen);
    ASSERT_EQ(tree.bv.getBit(5 * blockSize), closeParen);
    ASSERT_EQ(tree.bv.sizeInBits(), blockSize * 10);
    ASSERT_EQ(tree.numLeaves, 10);
    for (Index i = 0; i < 5; ++i) {
        ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(i)], 1 + i * blockSize);
        ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(i + 5)], 5 * blockSize - (i + 1) * blockSize) << i;
        if (i < 4) {
            ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(i) / 2], ((i / 2) * 2) * blockSize + 1) << i;
            ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(9 - i) / 2], ((i / 2) * 2) * blockSize) << i;
        }
    }
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(5) / 2], 4 * blockSize);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(0) / 4], 1);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(4) / 4], 2 * blockSize);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(9) / 4], 0);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(7) / 8], 1);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(9) / 8], 0);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(9) / 16], 0);
    ASSERT_EQ(tree.bitvecRmq(blockSize, 5 * blockSize), blockSize);
    ASSERT_EQ(tree.bitvecRmq(blockSize, str.size()), str.size() - 1);
    ASSERT_EQ(tree.bitvecRmq(0, blockSize), 0);
    ASSERT_EQ(tree.bitvecRmq(5 * blockSize - 6, 5 * blockSize + 4), 5 * blockSize - 6);
    ASSERT_EQ(tree.bitvecRmq(5 * blockSize - 4, 5 * blockSize + 6), 5 * blockSize + 5);
    testRmmTree(tree);
}

TEST(RmmTree, 9BlocksAscending) {
    Index blockSize = RangeMinMaxTree<>::blockSize;
    Index mid = 4 * blockSize + blockSize / 2;
    assert(blockSize % 2 == 0);
    std::string str(blockSize * 9, '1');
    for (Index i = str.size() / 2; i < str.size(); ++i) {
        str[i] = '0';
    }
    Allocation allocation(Bitvec::allocatedSizeInLimbsForBits(str.size()));
    RangeMinMaxTree<> tree((Bitvec(str)), allocation.memory());
    ASSERT_EQ(tree.bv.getBit(5 * blockSize), closeParen);
    ASSERT_EQ(tree.bv.getBit(4 * blockSize), openParen);
    ASSERT_EQ(tree.bv.getBit(mid), closeParen);
    ASSERT_EQ(tree.bv.getBit(mid - 1), openParen);
    ASSERT_EQ(tree.bv.sizeInBits(), blockSize * 9);
    ASSERT_EQ(tree.numLeaves, 9);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(4)], 4 * blockSize);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(5)], 3 * blockSize);
    for (Index i = 0; i < 4; ++i) {
        ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(i)], 1 + i * blockSize) << i;
        ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(8 - i)], i * blockSize) << i;
        ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(i) / 2], ((i / 2) * 2) * blockSize + 1) << i;
        if (i == 0) {
            ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(8) / 2], 0);
        } else {
            ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(8 - i) / 2], (((i - 1) / 2) * 2 + 1) * blockSize) << i;
        }
    }
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(4) / 2], 3 * blockSize);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(5) / 2], 3 * blockSize);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(0) / 4], 1);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(4) / 4], blockSize);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(8) / 4], 0);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(7) / 8], 1);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(8) / 8], 0);
    ASSERT_EQ(tree.rmmArr[tree.leafIdxInArr(9) / 16], 0);
    ASSERT_EQ(tree.bitvecRmq(blockSize, 5 * blockSize), blockSize);
    ASSERT_EQ(tree.bitvecRmq(blockSize, str.size()), str.size() - 1);
    ASSERT_EQ(tree.bitvecRmq(blockSize * 4, blockSize * 5), blockSize * 5 - 1);
    ASSERT_EQ(tree.bitvecRmq(0, blockSize), 0);
    ASSERT_EQ(tree.bitvecRmq(mid - 6, mid + 4), mid - 6);
    ASSERT_EQ(tree.bitvecRmq(mid - 4, mid + 6), mid + 5);
    testRmmTree(tree);
}

TEST(RmmTree, 10BlocksAlternating) {
    const Index blockSize = RangeMinMaxTree<>::blockSize;
    std::string str(blockSize * 10, '1');
    for (Index i = 0; i < blockSize; ++i) {
        str[i + 2 * blockSize] = '0';
        str[i + 4 * blockSize] = '0';
        str[i + 6 * blockSize] = '0';
        str[i + 8 * blockSize] = '0';
        str[i + 9 * blockSize] = '0';
    }
    Allocation allocation(Bitvec::allocatedSizeInLimbsForBits(str.size()));
    RangeMinMaxTree<> tree((Bitvec(str)), allocation.memory());
    ASSERT_EQ(tree.rmqImpl(blockSize, 5 * blockSize).minExcess, blockSize);
    ASSERT_EQ(tree.rmqImpl(blockSize - 1, blockSize + 3).minExcess, blockSize);
    ASSERT_EQ(tree.rmqImpl(2 * blockSize - 1, 2 * blockSize + 3).pos, 2 * blockSize + 2);
    ASSERT_EQ(tree.rmqImpl(2 * blockSize - 1, 2 * blockSize + 3).minExcess, 2 * blockSize - 3);
    assert(blockSize > 7);
    ASSERT_EQ(tree.rmqImpl(7 * blockSize + 2, 7 * blockSize + 7).minExcess, blockSize + 3);
    ASSERT_EQ(tree.bitvecRmq(blockSize, 2 * blockSize), blockSize);
    ASSERT_EQ(tree.bitvecRmq(2 * blockSize, 3 * blockSize), 3 * blockSize - 1);
    ASSERT_EQ(tree.bitvecRmq(blockSize, str.size()), str.size() - 1);
    ASSERT_EQ(tree.bitvecRmq(0, blockSize), 0);
    ASSERT_EQ(tree.bitvecRmq(7 * blockSize - 6, 7 * blockSize + 4), 7 * blockSize - 1);
    ASSERT_EQ(tree.bitvecRmq(7 * blockSize - 4, 7 * blockSize + 6), 7 * blockSize - 1);
    ASSERT_EQ(tree.bitvecRmq(8 * blockSize - 6, 8 * blockSize + 4), 8 * blockSize - 6);
    ASSERT_EQ(tree.bitvecRmq(8 * blockSize - 4, 8 * blockSize + 6), 8 * blockSize + 5);
    testRmmTree(tree);
}

TEST(RmmTree, IncreasingDegree) {
    std::string s("1");
    for (Index i = 0; i < 100; ++i) {
        for (Index j = 0; j < i; ++j) {
            s.push_back(('1'));
        }
        for (Index j = 0; j < i; ++j) {
            s.push_back(('0'));
        }
    }
    s.push_back('0');
    Allocation allocation(Bitvec::allocatedSizeInLimbsForBits(s.size()));
    RangeMinMaxTree<> tree((Bitvec(s)), allocation.memory());
    ASSERT_EQ(tree.rmqImpl(0, 2).minExcess, 1);
    ASSERT_EQ(tree.rmqImpl(0, 5000).minExcess, 1);
    ASSERT_EQ(tree.rmqImpl(1, 2345).minExcess, 1);
    ASSERT_EQ(tree.rmqImpl(100, 5432).minExcess, 1);
    ASSERT_EQ(tree.rmqImpl(1 + 10 * 11, 1 + 11 * 12).minExcess, 1);
    for (Index i = 0; i < 10; ++i) {
        ASSERT_EQ(tree.rmqImpl(1 + 10 * 11 + i, 1 + 11 * 12 - i).minExcess, 1 + i) << i;
        ASSERT_EQ(tree.rmqImpl(10 * 11, 1 + 11 * 12 - i).minExcess, 1) << i;
        ASSERT_EQ(tree.rmqImpl(1 + 10 * 11 + i, 1 + 11 * 12).minExcess, 1) << i;
    }
    testRmmTree(tree);
}


TEST(RmmTree, DecreasingDegreeAlternatingDegree1) {
    std::string s("11");
    std::vector<Index> startIndices;
    std::vector<Index> lengths;
    for (Index i = 100; i + 1 > 0; --i) {
        s.push_back('1');
        s.push_back('0');
        startIndices.push_back(s.size());
        lengths.push_back(i);
        for (Index j = 0; j < i; ++j) {
            s.push_back(('1'));
        }
        for (Index j = 0; j < i; ++j) {
            s.push_back(('0'));
        }
    }
    s.push_back('0');
    s.push_back('0');
    Allocation allocation(Bitvec::allocatedSizeInLimbsForBits(s.size()));
    RangeMinMaxTree<> tree((Bitvec(s)), allocation.memory());
    ASSERT_EQ(tree.rmqImpl(0, 2).minExcess, 1);
    ASSERT_EQ(tree.rmqImpl(0, 5000).minExcess, 1);
    ASSERT_EQ(tree.rmqImpl(1, 3).minExcess, 2);
    ASSERT_EQ(tree.rmqImpl(1, 2345).minExcess, 2);
    ASSERT_EQ(tree.rmqImpl(100, 5432).minExcess, 2);
    ASSERT_EQ(tree.rmqImpl(startIndices[10], startIndices[11]).minExcess, 2);
    for (Index i = 0; i < 20; ++i) {
        ASSERT_EQ(tree.rmqImpl(startIndices[20] + i, startIndices[21] - 2 - i).minExcess, 2 + i) << i;
        ASSERT_EQ(tree.rmqImpl(startIndices[20], startIndices[21] - 2 - i).minExcess, i == 0 ? 2 : 3) << i;
        ASSERT_EQ(tree.rmqImpl(startIndices[20] + i, startIndices[21] - 2).minExcess, 2) << i;
    }
    testRmmTree(tree);
}

TEST(RmmTree, Random) {
    SuccinctRMQ _ = randomRmq(12345);
    const RangeMinMaxTree<>& tree = _.getTree();
    const Index blockSize = RangeMinMaxTree<>::blockSize;
    ASSERT_EQ(tree.numLeaves, roundUpDiv(tree.bv.sizeInBits(), blockSize));
    ASSERT_EQ(tree.getBitvector().getBit(0), openParen);
    ASSERT_EQ(tree.getBitvector().getBit(tree.getBitvector().sizeInBits() - 1), closeParen);
    testRmmTree(tree);
    for (Index i = 0; i < tree.numLeaves; ++i) {
        ASSERT_EQ(tree.findMinInBlock(i * blockSize, (i + 1) * blockSize).minExcess, tree.rmmArr[tree.leafIdxInArr(i)]);
    }
}

#ifdef ADS_HAS_CPP20

TEST(RmmTree, Constexpr) {
    // TODO: Uncomment and implement
    //    static_assert(RangeMinMaxTree<>(Bitvec("10111000")).findMinInBlock(0, 3).pos == 1);
    //    static_assert(RangeMinMaxTree<>(Bitvec("10111000")).bitvecRmq(0, 8) == 7);
}

#endif
