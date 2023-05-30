#ifndef BITVECTOR_BITVEC_LAYOUTS_HPP
#define BITVECTOR_BITVEC_LAYOUTS_HPP

#include "common.hpp"
#include <cstring>
#include <iostream>

namespace ads {

#ifdef ADS_HAS_CPP20
template<typename T>
concept IsLayout = requires(T& t) {
    T();
    T(Index);
    { T::superblockSize() } -> std::convertible_to<Index>;
    { T::blockSize() } -> std::convertible_to<Index>;
    { T::numBlocksInSuperblock() } -> std::convertible_to<Index>;
    { T::bytesPerBlockCount() } -> std::convertible_to<Index>;
    { t.completeSizeInElems() } -> std::convertible_to<Index>;
    { t.getElem(Index()) } -> std::convertible_to<Elem&>;
    { t.getSuperblockCount(Index()) } -> std::convertible_to<Elem&>;
    { t.getBlockCount(Index(), Index()) } -> std::convertible_to<Index>;
    { t.setBlockCount(Index(), Index(), Index()) } -> std::same_as<void>;
    { t.numElems() } -> std::convertible_to<Index>;
    { t.numBlocks() } -> std::convertible_to<Index>;
    { t.numSuperblocks() } -> std::convertible_to<Index>;
    { T::completeSizeInElems(Index()) } -> std::convertibele_to<Index>;
};
#define ADS_LAYOUT_CONCEPT IsLayout
#else
#define ADS_LAYOUT_CONCEPT class
#endif

constexpr static Index dynSize = Index(-1);


template<Index NumBits>
struct BitwiseAccess {

    std::unique_ptr<Elem[]> vec;
    constexpr static Elem mask = (Elem(1) << NumBits) - 1;

    [[nodiscard]] Elem getBits(Index elemStartIndex, Index inElem) const noexcept {
        assert(elemStartIndex >= 0 && inElem >= 0 && NumBits > 0 && NumBits < 64);
        return (vec[elemStartIndex] >> (NumBits * inElem)) & mask;
    }

    void setBits(Index elemStartIndex, Index inElem, Elem value) noexcept {
        assert(elemStartIndex >= 0 && inElem >= 0 && NumBits > 0 && NumBits < 64);
        assert(value < (Elem(1) << NumBits));
        vec[elemStartIndex] &= ~(mask << (NumBits * inElem));
        vec[elemStartIndex] |= value << (NumBits * inElem);
    }
};

template<>
struct BitwiseAccess<dynSize> {

    std::unique_ptr<Elem[]> vec;

    [[nodiscard]] Elem getBits(Index elemStartIndex, Index bitStartIndex, Index numBits) const noexcept {
        assert(elemStartIndex >= 0 && bitStartIndex >= 0 && bitStartIndex < 64 && numBits > 0 && numBits < 64);
        Elem mask = (Elem(1) << numBits) - 1;
        Elem r = (vec[elemStartIndex] >> bitStartIndex) & mask;
        if (numBits + bitStartIndex > 64) {
            Index remaining = numBits + bitStartIndex - 64;
            mask = (Elem(1) << remaining) - 1;
            r += (vec[elemStartIndex + 1] & mask) << (numBits - remaining);
        }
        return r;
    }

    [[nodiscard]] Elem getBits(Index bitStartIndex, Index numBits) const noexcept {
        return getBits(bitStartIndex / 64, bitStartIndex % 64, numBits);
    }

    void setBitsImpl(Index elemStartIndex, Index bitStartIndex, Elem value, Index numBits) noexcept {
        assert(numBits > 0 && numBits < 64 && bitStartIndex >= 0 && bitStartIndex < 64);
        const Elem mask = (Elem(1) << numBits) - 1;
        vec[elemStartIndex] &= ~(mask << bitStartIndex);
        vec[elemStartIndex] |= (value & mask) << bitStartIndex;
    }

    void setBits(Index elemStartIndex, Index bitStartIndex, Elem value, Index numBits) noexcept {
        assert(elemStartIndex >= 0 && bitStartIndex >= 0 && numBits > 0 && numBits < 64);
        setBitsImpl(elemStartIndex, bitStartIndex, value, numBits);
        if (numBits + bitStartIndex > 64) {
            Index remaining = numBits + bitStartIndex - 64;
            setBitsImpl(elemStartIndex + 1, 0, value >> (numBits - remaining), remaining);
        }
    }

    void setBits(Index bitStartIndex, Elem value, Index numBits) noexcept {
        setBits(bitStartIndex / 64, bitStartIndex % 64, value, numBits);
    }
};

template<ADS_LAYOUT_CONCEPT Layout, Index NumBits>
struct DefaultLayoutImpl : BitwiseAccess<NumBits> {

    [[nodiscard]] constexpr const Layout& derived() const noexcept {
        return static_cast<const Layout&>(*this);
    }

    [[nodiscard]] constexpr static Index numBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, Layout::blockSize());
    }

    [[nodiscard]] constexpr Index numBlocks() const noexcept {
        return roundUpDiv(derived().numElems(), derived().blockSize());
    }

    [[nodiscard]] constexpr static Index numSuperblocks(Index numElems) noexcept {
        return roundUpDiv(numElems, Layout::superblockSize());
    }

    [[nodiscard]] constexpr Index numSuperblocks() const noexcept {
        return roundUpDiv(derived().numElems(), derived().superblockSize());
    }

    [[nodiscard]] constexpr static Index numBlocksInSuperblock() noexcept {
        assert(Layout::superblockSize() % Layout::blockSize() == 0);
        return Layout::superblockSize() / Layout::blockSize();
    }
};

// This layout wastes a lot of memory; it needs roughly double the total memory of the actual bit vector
struct CacheEfficientLayout : DefaultLayoutImpl<CacheEfficientLayout, 8> {

    using Base = DefaultLayoutImpl<CacheEfficientLayout, 8>;
    Index numElements;

    CacheEfficientLayout() noexcept : numElements(0) {}

    explicit CacheEfficientLayout(Index numElems) noexcept {
        vec = makeUniqueForOverwrite<Elem>(completeSizeInElems(numElems));
        numElements = numElems;
        // TODO: ensure that allocated memory is cache aligned
    }

    constexpr static Index superblockSize() noexcept {
        return 4;
    }

    constexpr static Index blockSize() noexcept {
        return 1;
    }

    constexpr static Index bytesPerBlockCount() noexcept {
        return 1;
    }

    static Index completeSizeInElems(Index numElements) noexcept {
        return roundUpDiv(numElements, 4) * ELEMS_PER_CACHELINE;
    }

    [[nodiscard]] Index numElems() const noexcept {
        return numElements;
    }

    static Index getElemIdx(Index i) noexcept {
        assert(i >= 0);
        return (i / superblockSize()) * (ELEMS_PER_CACHELINE) + i % superblockSize();// TODO: Make sure this gets optimized to bitshifts (possibly negative i may interfere)
    }

    Elem& getElem(Index i) noexcept {
        assert(i >= 0 && i < numElements);
        return vec[getElemIdx(i)];
    }
    [[nodiscard]] const Elem& getElem(Index i) const noexcept {
        assert(i >= 0 && i < numElements);
        return vec[getElemIdx(i)];
    }

    static Index getSuperblockCountIdx(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0);
        return superblockIdx * ELEMS_PER_CACHELINE + superblockSize();
    }

    [[nodiscard]] Elem& getSuperblockCount(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks());
        return vec[getSuperblockCountIdx(superblockIdx)];
    }
    [[nodiscard]] const Elem& getSuperblockCount(Index superblockIdx) const noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks());
        return vec[getSuperblockCountIdx(superblockIdx)];
    }

    Index getBlockCount(Index superblockIdx, Index blockIdx) const noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks() && blockIdx >= 0 && blockIdx < numBlocks());
        return Index((vec[getSuperblockCountIdx(superblockIdx) + 1] >> (8 * blockIdx)) & 0xff);
    }

    void setBlockCount(Index superblockIdx, Index blockIdx, std::uint8_t value) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks() && blockIdx >= 0 && blockIdx < numBlocks());
        Index i = getSuperblockCountIdx(superblockIdx) + 1;
        vec[i] &= ~(Elem(0xff) << (8 * blockIdx));
        vec[i] |= Elem(value) << (8 * blockIdx);
    }
};


/// TODO: Doxygen
template<Index SuperblockSize = (1 << 16) / 64, Index BlockSize = 1>
struct SimpleLayout : DefaultLayoutImpl<SimpleLayout<SuperblockSize, BlockSize>, BlockSize * 8> {

    static_assert(0 < BlockSize && BlockSize <= SuperblockSize);

private:
    ADS_CONSTEVAL static Index numBytesOfBlockCount(Index blocksInSuperblock) noexcept {
        // no constexpr std::bit_floor in C++17; using <= instead of < is fine because no entry actually stores this number
        return blocksInSuperblock <= 256 ? 1 : 1 + numBytesOfBlockCount(blocksInSuperblock / 256);
    }

public:
    using Base = DefaultLayoutImpl<SimpleLayout<SuperblockSize, BlockSize>, BlockSize * 8>;
    std::unique_ptr<Elem[]> vec = nullptr;
    Elem* ADS_RESTRICT superblocks = nullptr;
    unsigned char* ADS_RESTRICT blocks = nullptr;// examining the object representation is only valid through unsigned char* or std::byte*

    using Base::numBlocks;
    using Base::numBlocksInSuperblock;
    using Base::numSuperblocks;

    SimpleLayout() noexcept = default;

    explicit SimpleLayout(Index numElements) : vec(makeUniqueForOverwrite<Elem>(completeSizeInElems(numElements))) {
        superblocks = vec.get() + numElements;
        blocks = reinterpret_cast<unsigned char*>(superblocks + numSuperblocks());
    }

    [[nodiscard]] constexpr static Index superblockSize() noexcept {
        return SuperblockSize;
    }

    [[nodiscard]] constexpr static Index blockSize() noexcept {
        return BlockSize;
    }

    [[nodiscard]] constexpr static Index bytesPerBlockCount() noexcept {
        constexpr Index res = numBytesOfBlockCount(numBlocksInSuperblock());// force compile time evaluation in C++17 mode
        return res;
    }

    [[nodiscard]] Index completeSizeInElems() const noexcept {
        return completeSizeInElems(numElems());
    }

    [[nodiscard]] static Index completeSizeInElems(Index numElements) noexcept {
        return numElements + blockSize() * numBlocks(numElements) + superblockSize() * numSuperblocks(numElements);
    }

    [[nodiscard]] const Elem& getElem(Index i) const noexcept {
        return vec[i];
    }

    [[nodiscard]] Elem& getElem(Index i) noexcept {
        return vec[i];
    }

    [[nodiscard]] Elem& getSuperblockCount(Index i) const noexcept {
        return superblocks[i];
    }

    [[nodiscard]] Index getBlockCount(Index superblockIdx, Index blockIdx) const noexcept {
        constexpr static Index blockCountWidth = bytesPerBlockCount();
        unsigned char* ADS_RESTRICT ptr = blocks + (superblockIdx * numBlocksInSuperblock() + blockIdx) * blockCountWidth;
        if constexpr (blockCountWidth == 1) {
            return Index(*ptr);
        } else if constexpr (blockCountWidth == 2 || blockCountWidth == 4) {
            return ptrBitCast<IntType<blockCountWidth>>(ptr);
        } else if constexpr (blockCountWidth == 3) {
            static_assert(blockCountWidth == 3);
            return Elem(*ptr) + (Elem(ptr[1]) << 8) + (Elem(ptr[2]) << 16);// don't assume alignment
        }
    }

    void setBlockCount(Index superblockIdx, Index blockIdx, Index value) noexcept {
        constexpr static Index blockCountWidth = bytesPerBlockCount();
        auto ptr = blocks + (superblockIdx * numBlocksInSuperblock() + blockIdx) * blockCountWidth;
        assert(value >= 0 && value < Elem(1) << (blockCountWidth * 8));
        if constexpr (blockCountWidth == 1) {
            *ptr = value;
        } else if constexpr (blockCountWidth == 2 || blockCountWidth == 4) {
            IntType<blockCountWidth> tmp = value;
            std::memcpy(ptr, &value, sizeof(tmp));
        } else {
            static_assert(blockCountWidth == 3);
            *ptr = std::uint8_t(value);
            ptr[1] = std::uint8_t(value >> 8);
            ptr[2] = std::uint8_t(value >> 16);
        }
    }

    [[nodiscard]] Index numElems() const noexcept {
        return superblocks - vec.get();
    }
};

#ifdef ADS_HAS_CPP20
static_assert(Layout<CacheEfficientLayout>);
static_assert(Layout<NaiveLayout4ElemSuperblocks>);
static_assert(Layout<NaiveLayout1024ElemSuperblocks>);
#endif

}// namespace ads

#endif//BITVECTOR_BITVEC_LAYOUTS_HPP
