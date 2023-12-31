#include "../include/elias_fano.hpp"
#include "gtest/gtest.h"
#include <random>

using namespace ads;

TEST(EliasFano, ConstructionSmallAscending) {
    Elem arr[] = {1, 2, 3};
    EliasFano<> ef(arr);
    ASSERT_EQ(ef.size(), 3);
    ASSERT_EQ(ef.numLowerBitsPerNumber(), 0);
    ASSERT_EQ(ef.numBitsPerNumber(), 3);
    const auto& upper = ef.getUpper();
    ASSERT_LE(upper.size(), 12);
    ASSERT_TRUE(upper.getBit(0));
    for (Index i = 0; i < 3; ++i) {
        ASSERT_FALSE(upper.getBit(2 * i + 1)) << i;
    }
    ASSERT_EQ(upper.selectOne(0), 0);
    ASSERT_EQ(upper.selectOne(1), 2);
    ASSERT_EQ(upper.selectZero(0), 1);
    ASSERT_EQ(upper.selectZero(2), 5);
    ASSERT_EQ(ef.getSmallest(), 1);
    //    ASSERT_EQ(upper, EfficientBitvec<>("010101"));
    ASSERT_EQ(upper, EfficientBitvec<>("1010101"));
    const auto& lower = ef.getLower();
    ASSERT_EQ(lower.bitAccess.numBits, 0);
    for (Index i = 0; i < std::size(arr); ++i) {
        ASSERT_EQ(ef.get(i), i + 1) << i;
    }
}

TEST(EliasFano, ConstructionSmall) {
    EliasFano<> ef{12, 14, 990, 991, 100'000};
    ASSERT_EQ(ef.size(), 5);
    ASSERT_GT(ef.numLowerBitsPerNumber(), 0);
    ASSERT_LT(ef.numLowerBitsPerNumber(), std::ceil(intLog2(100'000u)));
    ASSERT_EQ(ef.getLower()[0] % 2, 0);
    ASSERT_EQ(ef.getLower()[3] % 2, 1);
    ASSERT_EQ(ef.getLower()[0], 0);
    ASSERT_EQ(ef.getLower()[1], 2 & ((1 << ef.numLowerBitsPerNumber()) - 1));
    ASSERT_EQ(ef.getLower()[2], 978 & ((1 << ef.numLowerBitsPerNumber()) - 1));
    ASSERT_EQ(ef.get(0), 12);
    ASSERT_EQ(ef.get(1), 14);
    ASSERT_EQ(ef.get(2), 990);
    ASSERT_EQ(ef.get(3), 991);
    ASSERT_EQ(ef.get(4), 100'000);
    ASSERT_EQ(ef.getLower().bitAccess.numBits, ef.numLowerBitsPerNumber());
}

TEST(EliasFano, ConstructionSmallNegative) {
    EliasFano<Index> ef{-100, -5, 4, 12};
    ASSERT_EQ(ef.size(), 4);
    ASSERT_EQ(ef.get(0), -100);
    ASSERT_EQ(ef.get(1), -5);
    ASSERT_EQ(ef.get(2), 4);
    ASSERT_EQ(ef.get(3), 12);
    ASSERT_LE(ef.numLowerBitsPerNumber(), 7);
}

TEST(EliasFano, Empty) {
    EliasFano<> ef{};
    ASSERT_EQ(ef.size(), 0);
}

TEST(EliasFano, OneElement) {
    EliasFano<> ef{42};
    ASSERT_EQ(ef.size(), 1);
    ASSERT_EQ(ef.get(0), 42);
    ASSERT_EQ(ef.predecessor(44), 42);
    ASSERT_LE(ef.numBitsPerNumber(), 2);
    ASSERT_LE(ef.numUpperBitsPerNumber(), 2);
}

void testElems(const std::vector<Elem>& arr, const EliasFano<>& ef, std::string_view name) {
    assert(std::is_sorted(arr.begin(), arr.end()));
    ASSERT_EQ(ef.size(), ssize(arr)) << name;
    Index numLowerBits = ef.numLowerBitsPerNumber();
    ASSERT_LE(numLowerBits, 64 - std::ceil(std::log2(ef.size() + 2)));
    ASSERT_GE(numLowerBits, 0);
    ASSERT_EQ(ef.getSmallest(), arr[0]);
    ASSERT_TRUE(ef.getUpper().getBit(ef.getUpper().size() - 1));
    for (Index i = 0; i < ssize(arr); ++i) {
        if (ef.numLowerBitsPerNumber() > 0) {
            ASSERT_EQ(ef.getLower()[i], (arr[i] - arr[0]) % (Elem(1) << numLowerBits))
                    << i << " of " << ssize(arr) << ", name " << name;
        }
        ASSERT_EQ(ef[i], arr[i]) << i << " of " << ssize(arr) << ", name " << name;
    }
    ASSERT_TRUE(std::equal(arr.begin(), arr.end(), ef.numbers().begin()));
}

TEST(EliasFano, ConstructionLarge) {
    std::vector<Elem> arr((1 << 12) - 1);
    std::iota(arr.begin(), arr.end(), 0);
    EliasFano<> ef(arr);
    testElems(arr, ef, "dense");
    for (Index i = 0; i < ssize(arr); i += 2) {
        arr[i] += i * i + 1;
    }
    arr.push_back(0);
    std::sort(arr.begin(), arr.end());
    testElems(arr, EliasFano<>(arr), "duplicate elements, power of two");
    arr.push_back(Elem(-1));
    assert(std::is_sorted(arr.begin(), arr.end()));
    ef = EliasFano<>(arr);
    ASSERT_EQ(ef.getLower()[ssize(arr) - 1], (Elem(1) << ef.numLowerBitsPerNumber()) - 1);
    testElems(arr, ef, "power of two + 1");
}

TEST(EliasFano, ConstructionRandom) {
    std::mt19937_64 engine(createRandomEngine());
    std::uniform_int_distribution<Elem> dist;
    std::vector<Elem> elems(dist(engine) & ((1 << 16) - 1));
    for (auto& elem : elems) {
        elem = dist(engine);
    }
    std::sort(elems.begin(), elems.end());
    EliasFano ef(elems);
    testElems(elems, ef, "random");
}

TEST(EliasFano, Select) {
    EliasFano<> ef{0, 3, 3, Elem(-1)};
    ASSERT_EQ(ef.numBitsPerNumber(), 64);
    ASSERT_EQ(ef.numBitsPerNumber() - ef.numLowerBitsPerNumber(), 3);
    //    ASSERT_EQ(ef.getUpper(), EfficientBitvec<>("011100000001"));
    ASSERT_EQ(ef.getUpper(), EfficientBitvec<>("1000111111101"));
    ASSERT_EQ(ef.get(0), 0);
    ASSERT_EQ(ef.get(1), 3);
    ASSERT_EQ(ef.get(2), 3);
    ASSERT_EQ(ef.get(3), Elem(-1));
}

TEST(EliasFano, PredecessorSmall) {
    EliasFano<> ef{1, 4, 5, 6, 10, Elem(-3)};
    ASSERT_EQ(ef.numBitsPerNumber(), 64);
    ASSERT_EQ(ef.numBitsPerNumber() - ef.numLowerBitsPerNumber(), 3);
    ASSERT_EQ(ef.predecessor(1), 1);
    ASSERT_EQ(ef.predecessor(2), 1);
    ASSERT_EQ(ef.predecessor(3), 1);
    ASSERT_EQ(ef.predecessor(4), 4);
    ASSERT_EQ(ef.predecessor(5), 5);
    for (Index i = 6; i < 10; ++i) {
        ASSERT_EQ(ef.predecessor(i), 6);
    }
    for (Index i = 10; i < 10'000; ++i) {
        ASSERT_EQ(ef.predecessor(i), 10);
    }
    ASSERT_EQ(ef.predecessor(Elem(-1)), Elem(-3));
}

TEST(EliasFano, GetSmallSizeBigValues) {
    EliasFano<> ef{Elem(-10), Elem(-3), Elem(-3), Elem(-1)};
    ASSERT_EQ(ef.get(0), Elem(-10));
    ASSERT_EQ(ef.get(1), Elem(-3));
    ASSERT_EQ(ef.get(2), ef.get(1));
    ASSERT_EQ(ef.get(3), Elem(-1));
    ASSERT_EQ(ef.predecessor(Elem(-10)), Elem(-10));
    ASSERT_EQ(ef.predecessor(Elem(-4)), Elem(-10));
    ASSERT_EQ(ef.predecessor(Elem(-3)), Elem(-3));
    ASSERT_EQ(ef.predecessor(Elem(-2)), Elem(-3));
    ASSERT_EQ(ef.predecessor(Elem(-1)), Elem(-1));
}

TEST(EliasFano, Small) {
    EliasFano<> ef{1234};
    ASSERT_EQ(ef.get(0), 1234);
    ASSERT_EQ(ef.predecessor(123456), 1234);
    ef = EliasFano<>{0};
    ASSERT_EQ(ef.get(0), 0);
    ASSERT_EQ(ef.predecessor(0), 0);
    ASSERT_EQ(ef.predecessor(1), 0);
    ef = EliasFano<>{Elem(-2)};
    ASSERT_EQ(ef.get(0), Elem(-2));
    ASSERT_EQ(ef.predecessor(Elem(-1)), Elem(-2));
    ef = EliasFano<>{1, 2};
    ASSERT_EQ(ef.get(0), 1);
    ASSERT_EQ(ef.get(1), 2);
    ASSERT_EQ(ef.predecessor(2), 2);
    ASSERT_EQ(ef.predecessor(10), 2);
    ASSERT_EQ(ef.predecessor(1), 1);
    ef = EliasFano<>{4, Elem(-3)};
    ASSERT_EQ(ef.get(0), 4);
    ASSERT_EQ(ef.get(1), Elem(-3));
    ASSERT_EQ(ef.predecessor(4), 4);
    ASSERT_EQ(ef.predecessor(Elem(-4)), 4);
    ASSERT_EQ(ef.predecessor(Elem(-1)), Elem(-3));
    ef = EliasFano<>{1, 1};
    ASSERT_EQ(ef.get(0), 1);
    ASSERT_EQ(ef.get(1), 1);
    ASSERT_EQ(ef.predecessor(1), 1);
    ASSERT_EQ(ef.predecessor(1234567), 1);
    ef = EliasFano<>{9, 9, Elem(-1) / 2};
    ASSERT_EQ(ef.get(0), 9);
    ASSERT_EQ(ef.get(1), 9);
    ASSERT_EQ(ef.get(2), Elem(-1) / 2);
    ASSERT_EQ(ef.predecessor(9), 9);
    ASSERT_EQ(ef.predecessor(10), 9);
    ASSERT_EQ(ef.predecessor(Elem(-1) / 2), Elem(-1) / 2);
    ASSERT_EQ(ef.predecessor(Elem(-1) / 2 + 1), Elem(-1) / 2);
    ASSERT_EQ(ef.predecessor(Elem(-1) - 1), Elem(-1) / 2);
    ef = EliasFano<>{Elem(-4), Elem(-2)};
    ASSERT_EQ(ef.get(0), Elem(-4));
    ASSERT_EQ(ef.get(1), Elem(-2));
    ASSERT_EQ(ef.predecessor(Elem(-1)), Elem(-2));
    ASSERT_EQ(ef.successor(1), Elem(-4));
}

TEST(EliasFano, SmallSignedChar) {
    EliasFano<signed char> charEF{-12, -3, 8, 13, 99};
    ASSERT_EQ(charEF.size(), 5);
    ASSERT_EQ(charEF.get(0), -12);
    ASSERT_EQ(charEF.get(3), 13);
    ASSERT_EQ(charEF.get(2), 8);
    ASSERT_EQ(charEF.predecessor(127), 99);
    ASSERT_EQ(charEF.successor(-2), 8);
    ASSERT_EQ(charEF.successor(12), 13);
    ASSERT_EQ(charEF.predecessor(-4), -12);
    charEF = EliasFano<signed char>{-128};
    ASSERT_EQ(charEF.size(), 1);
    ASSERT_EQ(charEF.get(0), -128);
    ASSERT_EQ(charEF.predecessor(127), -128);
    charEF = EliasFano<signed char>{-128, 127};
    ASSERT_EQ(charEF.size(), 2);
    ASSERT_EQ(charEF.get(1), 127);
    ASSERT_EQ(charEF.predecessor(123), -128);
    ASSERT_EQ(charEF.successor(-12), 127);
}

TEST(EliasFano, ManyIdentical) {
    std::vector<Elem> values(10'000, 1234567890);
    EliasFano<> ef(values);
    ASSERT_EQ(ef.size(), 10'000);
    ASSERT_EQ(ef.get(0), 1234567890);
    ASSERT_EQ(ef.predecessor(1234567891), 1234567890);
    ASSERT_EQ(ef.successor(0), 1234567890);
    values.insert(values.begin(), 0);
    ef = EliasFano<>(values);
    ASSERT_EQ(ef.size(), 10'001);
    ASSERT_EQ(ef.get(0), 0);
    ASSERT_EQ(ef.get(1), 1234567890);
    ASSERT_EQ(ef.successor(0), 0);
    ASSERT_EQ(ef.successor(1), 1234567890);
    values = std::vector<Elem>(12'345, 0);
    ef = EliasFano<>(values);
    ASSERT_EQ(ef.size(), 12345);
    ASSERT_EQ(ef.get(0), 0);
    ASSERT_EQ(ef.predecessor(Elem(-1)), 0);
    values.push_back(10);
    ef = EliasFano<>(values);
    ASSERT_EQ(ef.get(0), 0);
    ASSERT_EQ(ef.get(1), 0);
    ASSERT_EQ(ef.get(values.size() - 1), 10);
    ASSERT_EQ(ef.predecessor(1), 0);
    ASSERT_EQ(ef.predecessor(11), 10);
}


TEST(EliasFano, PredecessorAscending) {
    std::vector<Elem> vec(1'000'000);
    std::iota(vec.begin(), vec.end(), 0);
    EliasFano<> ef(vec);
    ASSERT_LE(ef.numBitsPerNumber(), intLog2(ssize(vec)) + 1);
    ASSERT_EQ(ef.numLowerBitsPerNumber(), 0);
    ASSERT_LE(ef.numAllocatedBits(), 2 * (ssize(vec) + ssize(vec) / 5));
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist(0, vec.back() + vec.back() / 100);
    for (Index i = 0; i < 10'000; ++i) {
        Elem searched = dist(engine);
        Elem pred = searched > vec.back() ? vec.back() : searched;
        ASSERT_EQ(pred, ef.predecessor(searched));
    }
}

TEST(EliasFano, PredecessorSparseNonexisting) {
    std::vector<Elem> vec(1'000'000);
    std::iota(vec.begin(), vec.end(), 65535);
    for (Index i = 0; i < ssize(vec); i += 1 / 100 + i % 4 + 1) {
        --vec[i];
    }
    EliasFano<> ef(vec);
    ASSERT_LE(ef.numBitsPerNumber(), intLog2(ssize(vec)) + 1);
    ASSERT_EQ(ef.numLowerBitsPerNumber(), 0);
    ASSERT_LE(ef.numAllocatedBits(), 2 * (ssize(vec) + ssize(vec) / 5));
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist(0, vec.back() + vec.back() / 100);
    for (Index i = 0; i < 10'000; ++i) {
        Elem searched = dist(engine) + ef.getSmallest();
        Elem expected = *(std::upper_bound(vec.begin(), vec.end(), searched) - 1);
        ASSERT_EQ(ef.predecessor(searched), expected);
    }
}

TEST(EliasFano, PredecessorRandom) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist;
    std::vector<Elem> vec;
    for (Index i = 0; i < 1'000'000; ++i) {
        vec.push_back(dist(engine));
    }
    std::sort(vec.begin(), vec.end());
    EliasFano<> ef(vec);
    for (Index i = 0; i < 100'000; ++i) {
        Elem searched = dist(engine);
        auto iter = std::upper_bound(vec.begin(), vec.end(), searched);
        if (iter == vec.begin()) {
            continue;
        }
        Elem pred = *(iter - 1);
        ASSERT_EQ(pred, ef.predecessor(searched));
    }
}

TEST(EliasFano, PredecessorRandomSmallIntervall) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist(0, 1'200'000);
    std::vector<Elem> vec;
    for (Index i = 0; i < 1'000'000; ++i) {
        vec.push_back(dist(engine));
    }
    std::sort(vec.begin(), vec.end());
    EliasFano<> ef(vec);
    for (Index i = 0; i < 100'000; ++i) {
        Elem searched = dist(engine);
        auto iter = std::upper_bound(vec.begin(), vec.end(), searched);
        if (iter == vec.begin()) {
            continue;
        }
        Elem pred = *(iter - 1);
        ASSERT_EQ(pred, ef.predecessor(searched));
    }
}

TEST(EliasFano, SuccessorRandom) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist;
    std::vector<Elem> vec;
    for (Index i = 0; i < 1'000'000; ++i) {
        vec.push_back(dist(engine));
    }
    std::sort(vec.begin(), vec.end());
    EliasFano<> ef(vec);
    ASSERT_EQ(ef.getLargest(), vec.back());
    for (Index i = 0; i < 10'000; ++i) {
        Elem searched = dist(engine);
        auto iter = std::lower_bound(vec.begin(), vec.end(), searched);
        if (iter == vec.end()) {
            continue;
        }
        Elem suc = *iter;
        ASSERT_EQ(suc, ef.successor(searched));
    }
}

TEST(EliasFano, Clustered) {
    std::vector<Limb> vec;
    for (Index i = 0; i < 1000; ++i) {
        vec.push_back(i);
    }
    for (Index i = 250'000; i < 251'000; ++i) {
        vec.push_back(i);
    }
    for (Index i = 500'000'000; i < 500'001'000; ++i) {
        vec.push_back(i);
    }
    for (Index i = 750'000'000'000; i < 750'000'001'000; ++i) {
        vec.push_back(i);
    }
    EliasFano<> ef(vec);
    for (auto number : vec) {
        ASSERT_EQ(ef.predecessor(number), number);
    }
}


#ifdef ADS_HAS_CPP20
// Constexpr tests implicitly test for UB -- the compiler is required to emit an error in that case
TEST(EliasFano, Constexpr) {
    using T = EliasFano<>;
    static_assert(T{1, 2, 3}.size() == 3);
    static_assert(T{4, 5, 6}.get(1) == 5);
    static_assert(T{0, 1234, 9999}.predecessor(9000) == 1234);
    static_assert(T{0, Elem(-1)}.successor(99999999) == Elem(-1));
    static_assert(T{}.size() == 0);
    static_assert(T{0}.get(0) == 0);
    static_assert(EliasFano<signed char>{-128, 127}.predecessor(0) == -128);
    static_assert(EliasFano<signed char>(std::vector<signed char>(42, 5)).predecessor(43) == 5);
}
#endif
