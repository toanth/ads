#include "../include/elias_fano.hpp"
#include "gtest/gtest.h"
#include <random>

using namespace ads;

TEST(EliasFano, ConstructionSmall) {
    Elem arr[] = {1, 2, 3};
    EliasFano<> ef(arr);
    ASSERT_EQ(ef.size(), 3);
    ASSERT_EQ(ef.numLowerBitsPerNumber(), 62);
    const auto& upper = ef.getUpper();
    ASSERT_EQ(upper.sizeInBits(), 3 * 2 + 2);
    ASSERT_EQ(upper.getBit(0), false);
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(upper.getBit(i + 1), true) << i;
    }
    ASSERT_EQ(upper.getBit(30), false);
    ASSERT_EQ(upper.selectZero(0), 0);
    ASSERT_EQ(upper.selectZero(1), 4);
    ASSERT_EQ(upper.selectOne(0), 1);
    ASSERT_EQ(upper.selectOne(2), 3);
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(ef.getUpperPart(i), 0) << i;
    }
    ASSERT_EQ(upper, Bitvector<>("01110000"));
    const auto& lower = ef.getLower();
    for (Index i = 0; i < 3; ++i) {
        ASSERT_EQ(lower.getBits(i), i + 1);
        ASSERT_EQ(ef.getLowerPart(i), i + 1);
    }
    for (Index i = 0; i < std::size(arr); ++i) {
        ASSERT_EQ(ef.get(i), i + 1) << i;
    }
}

TEST(EliasFano, Empty) {
    EliasFano<> ef{};
    ASSERT_EQ(ef.size(), 0);
}

void testElems(const std::vector<Elem>& arr, const EliasFano<>& ef, std::string_view name) {
    assert(std::is_sorted(arr.begin(), arr.end()));
    ASSERT_EQ(ef.size(), arr.size());
    Index numLowerBits = ef.numLowerBitsPerNumber();
    ASSERT_EQ(numLowerBits, 64 - std::ceil(std::log2(ef.size())));
    for (Index i = 0; i < arr.size(); ++i) {
        ASSERT_EQ(ef.getUpperPart(i), arr[i] >> numLowerBits) << i << " of " << arr.size() << ", name " << name;
        ASSERT_EQ(ef.getLowerPart(i), arr[i] % (Elem(1) << numLowerBits)) << i << " of " << arr.size() << ", name " << name;
        ASSERT_EQ(ef[i], arr[i]) << i << " of " << arr.size() << ", name " << name;
    }
    ASSERT_TRUE(std::equal(arr.begin(), arr.end(), ef.numbers().begin()));
}

TEST(EliasFano, ConstructionLarge) {
    std::vector<Elem> arr((1 << 12) - 1);
    std::iota(arr.begin(), arr.end(), 0);
    EliasFano<> ef(arr);
    testElems(arr, ef, "dense");
    for (Index i = 0; i < arr.size(); i += 2) {
        arr[i] += i * i + 1;
    }
    arr.push_back(0);
    std::sort(arr.begin(), arr.end());
    testElems(arr, EliasFano<>(arr), "duplicate elements");
    arr.push_back(Elem(-1));
    std::sort(arr.begin(), arr.end());
    ef = EliasFano<>(arr);
    ASSERT_EQ(ef.getUpperPart(arr.size() - 1), Elem(-1) >> ef.numLowerBitsPerNumber());
    ASSERT_EQ(ef.getLowerPart(arr.size() - 1), (Elem(1) << ef.numLowerBitsPerNumber()) - 1);
    testElems(arr, ef, "large-ish");
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
// TODO: Test with very few very large elements

TEST(EliasFano, Select) {
    EliasFano<> ef{0, 3, 3, Elem(-1)};
    ASSERT_EQ(ef.numUpperBitsPerNumber(), 2);
    ASSERT_EQ(ef.getUpper(), Bitvector<>("011100010"));
}

TEST(EliasFano, PredecessorSmall) {
    EliasFano<> ef{1, 4, 5, 6, 10, Elem(-3)};
    ASSERT_EQ(ef.numUpperBitsPerNumber(), 3);
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

TEST(EliasFano, PredecessorRandom) {
    auto engine = createRandomEngine();
    std::uniform_int_distribution<Elem> dist;
    std::vector<Elem> vec;
    for (Index i = 0; i < 1'000'000; ++i) {
        vec.push_back(dist(engine));
    }
    std::sort(vec.begin(), vec.end());
    EliasFano<> ef(vec);
    for (Index i = 0; i < 10'000; ++i) {
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
