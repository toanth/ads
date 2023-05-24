#ifndef BITVECTOR_BITVEC_LAYOUTS_HPP
#define BITVECTOR_BITVEC_LAYOUTS_HPP

#include "common.hpp"
#include <iostream>

namespace ads {

using Elem = std::uint64_t;
// This layout wastes a lot of memory; it needs roughly roughly double the memory in total than the actual bit vector
struct CacheEfficientLayout {

    CacheEfficientLayout([[maybe_unused]] Index numElems) {}

    constexpr static Index superblockSize() noexcept {
        return 4;
    }

    static Index completeSizeInElems(Index numElements) noexcept {
        return (numElements + 3) / 4 * ELEMS_PER_CACHELINE;
    }

    static Index getElemIdx(Index i) noexcept {
        assert(i >= 0);
        return (i / superblockSize()) * (ELEMS_PER_CACHELINE) + i % superblockSize();// TODO: Make sure this gets optimized to bitshifts (possibly negative i may interfere)
    }

    static Index getSuperblockCountIdx(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0);
        return superblockIdx * ELEMS_PER_CACHELINE + superblockSize();
    }

    static Index getBlockCount(Elem* ptr, Index superblockIdx, Index blockIdx) noexcept {
        assert(superblockIdx >= 0 && blockIdx >= 0);
        return Index((ptr[superblockIdx * ELEMS_PER_CACHELINE + superblockSize() + 1] >> (8 * blockIdx)) & 0xff);
    }

    static void setBlockCount(Elem* ptr, Index superblockIdx, Index blockIdx, std::uint8_t value) noexcept {
        assert(superblockIdx >= 0 && blockIdx >= 0);
        Index i = superblockIdx * ELEMS_PER_CACHELINE + superblockSize() + 1;
        ptr[i] &= ~(Elem(0xff) << (8 * blockIdx));
        ptr[i] |= Elem(value) << (8 * blockIdx);
    }
};

// 4 element superblocks make this similar to the cache efficient layout but use less space, roughly 3/8th additional memory
struct NaiveLayout4ElemSuperblocks {

    Index numElems;

    constexpr explicit NaiveLayout4ElemSuperblocks(Index numElems) noexcept : numElems(numElems) {}

    static Index numBlocks(Index numElements) noexcept {
        return numSuperblocks(numElements) * superblockSize();
    }

    constexpr static Index superblockSize() noexcept {
        //        Index s = blockSize(numElements);
        //        return s * s;
        return 4;// 4 Elems are 1 superblock
    }

    static Index numSuperblocks(Index numElems) noexcept {
        Index sPrime = superblockSize();
        return (numElems + sPrime - 1) / sPrime;
    }

    [[nodiscard]] Index numSuperblocks() const noexcept {
        return numSuperblocks(numElems);
    }

    static Index completeSizeInElems(Index numElements) noexcept {
        return numElements + (numBlocks(numElements) + 7) / 8 + numSuperblocks(numElements);
    }

    [[nodiscard]] Index getElemIdx(Index i) const noexcept {
        assert(i >= 0 && i <= numElems);
        return i;
    }

    [[nodiscard]] Index getSuperblockCountIdx(Index superblockIdx) const noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks(numElems));
        return numElems + superblockIdx;
    }

    Index getBlockCount(const Elem* ptr, Index superblockIdx, Index blockIdx) const noexcept {
        assert(superblockIdx >= 0 && blockIdx >= 0 && superblockIdx < numElems && blockIdx < superblockSize());
        Index blockToByte = superblockIdx * superblockSize() + blockIdx;
        Index arrIndex = numElems + numSuperblocks() + blockToByte / 8;
        assert(arrIndex >= 0 && arrIndex < completeSizeInElems(numElems));
        return Index((ptr[arrIndex] >> ((blockToByte % 8) * 8)) & 0xff);// TODO: Store blockBeginPtr? Measure!
    }

    void setBlockCount(Elem* ptr, Index superblockIdx, Index blockIdx, std::uint8_t value) const noexcept {
        Index i = superblockIdx * superblockSize() + blockIdx;
        Index j = numElems + (numElems + superblockSize() - 1) / superblockSize() + i / 8;
        Index shift = (i % 8) * 8;
        ptr[j] &= ~(Elem(0xff) << shift);
        ptr[j] |= Elem(value) << shift;
    }
};


struct NaiveLayout1024ElemSuperblocks {

    Index numElems;

    constexpr explicit NaiveLayout1024ElemSuperblocks(Index numElems) noexcept : numElems(numElems) {}

    static Index numBlocks(Index numElements) noexcept {
        return numSuperblocks(numElements) * superblockSize();
    }

    constexpr static Index superblockSize() noexcept {
        // 1024 Elems = 8192 Bytes = 65536 bits = 2^16 bits, which means 16bit variables can represent the block count
        return 1024;
    }

    static Index numSuperblocks(Index numElements) noexcept {
        Index sPrime = superblockSize();
        return (numElements + sPrime - 1) / sPrime;
    }

    [[nodiscard]] Index numSuperblocks() const noexcept {
        return numSuperblocks(numElems);
    }

    static Index completeSizeInElems(Index numElements) noexcept {
        return numElements + (numBlocks(numElements) + 3) / 4 + numSuperblocks(numElements);
    }

    [[nodiscard]] Index getElemIdx(Index i) const noexcept {
        assert(i >= 0 && i <= numElems);
        return i;
    }

    [[nodiscard]] Index getSuperblockCountIdx(Index superblockIdx) const noexcept {
        assert(superblockIdx >= 0);
        return numElems + superblockIdx;
    }

    Index getBlockCount(const Elem* ptr, Index superblockIdx, Index blockIdx) const noexcept {
        assert(superblockIdx >= 0 && blockIdx >= 0);
        Index blockNum = superblockIdx * superblockSize() + blockIdx;
        return Index((ptr[numElems + numSuperblocks() + blockNum / 4] >> ((blockNum % 4) * 16)) & 0xff'ff);// TODO: Store blockBeginPtr? Measure!
    }

    void setBlockCount(Elem* ptr, Index superblockIdx, Index blockIdx, std::uint16_t value) const noexcept {
        Index blockNum = superblockIdx * superblockSize() + blockIdx;
        Index i = numElems + numSuperblocks() + blockNum / 4;
        Index shift = (blockNum % 4) * 16;
        ptr[i] &= ~(Elem(0xff'ff) << shift);
        ptr[i] |= Elem(value) << shift;
    }
};

}// namespace ads

#endif//BITVECTOR_BITVEC_LAYOUTS_HPP
