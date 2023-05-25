#include "../include/bitvector.hpp"
#include <iostream>

int main() {
    std::cout << "Hello, World!" << std::endl;

    using namespace ads;
    Index numBits = 3;
    Index remainingBits = numBits;
    Bitvector<CacheEfficientLayout> bv(numBits);
    assert(bv.sizeInBits() == numBits);
    for (Index i = 0; i < bv.numSuperblocks(); ++i) {
        auto s = bv.superblockElems(i);
        assert(s.size() == bv.superblockSize());
        for (Index j = 0; j < bv.superblockSize(); ++j) {
            assert(s[j] == bv.element(i * bv.superblockSize() + j));
            if (remainingBits >= 64) {
                s[j] = Elem(-1);
                remainingBits -= 64;
            } else if (remainingBits == 0) {
                s[j] = 0;
            } else {
                s[j] = Elem(-1) << (64 - remainingBits);
                remainingBits = 0;
            }
        }
        bv.buildRankMetadata(i);
    }
    std::cout << "sizeInBits: " << bv.sizeInBits() << ", size in elems: " << bv.sizeInElems() << ", num superblocks: " << bv.numSuperblocks() << std::endl;
    for (Index i = 0; i < bv.sizeInElems(); ++i) {
        std::cout << std::hex << bv.element(i) << ' ';
    }
    std::cout << std::dec << std::endl;
    for (Index i = 0; i < bv.numSuperblocks(); ++i) {
        std::cout << "superblock " << i << ": " << bv.superblockCount(i) << '\n';
        for (Index j = 0; j < bv.superblockSize(); ++j) {
            std::cout << bv.blockCount(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "rank queries:" << std::endl;
    for (Index i = 0; i < numBits; ++i) {
        std::cout << bv.rankOne(i) << ' ';
    }

    return 0;
}
