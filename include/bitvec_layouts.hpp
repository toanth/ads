#ifndef BITVECTOR_BITVEC_LAYOUTS_HPP
#define BITVECTOR_BITVEC_LAYOUTS_HPP

#include "bit_access.hpp"
#include "common.hpp"
#include <cstring>
#include <iostream>

namespace ads {

#ifdef ADS_HAS_CPP20
template<typename T>
concept IsLayout = requires(T& t, const T& ct) {
    T();
    T(Index());
    { T::superblockSize() } -> std::convertible_to<Index>;
    { T::blockSize() } -> std::convertible_to<Index>;
    { T::numBlocksInSuperblock() } -> std::convertible_to<Index>;
    { T::bytesPerBlockCount() } -> std::convertible_to<Index>;
    { ct.allocatedSizeInElems() } -> std::convertible_to<Index>;
    { ct.getElem(Index()) } -> std::convertible_to<const Elem&>;
    { ct.getSuperblockCount(Index()) } -> std::convertible_to<const Elem&>;
    t.setSuperblockCount(Index(), Elem());
    { ct.getBlockCount(Index()) } -> std::convertible_to<Index>;
    t.setBlockCount(Index(), Index(), Index());
    { ct.numElems() } -> std::convertible_to<Index>;
    { ct.numBlocks() } -> std::convertible_to<Index>;
    { ct.numSuperblocks() } -> std::convertible_to<Index>;
    { T::allocatedSizeInElems(Index()) } -> std::convertible_to<Index>;
};
#define ADS_LAYOUT_CONCEPT IsLayout
#else
#define ADS_LAYOUT_CONCEPT class
#endif

template<typename Layout, Index NumBits>
struct DefaultLayoutImpl : BitStorage<NumBits> {

    [[nodiscard]] constexpr const Layout& derived() const noexcept { return static_cast<const Layout&>(*this); }

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

    [[nodiscard]] ADS_CONSTEVAL static Index numBlocksInSuperblock() noexcept {
        assert(Layout::superblockSize() % Layout::blockSize() == 0);
        return Layout::superblockSize() / Layout::blockSize();
    }

    [[nodiscard]] ADS_CONSTEVAL static Index bytesPerBlockCount() noexcept {
        return bytesNeededForIndexing(numBlocksInSuperblock());
    }
};

// This layout wastes a lot of memory; it needs roughly double the total memory of the actual bit vector
// Still, for small bitvectors, this may be a good option
struct CacheEfficientLayout : DefaultLayoutImpl<CacheEfficientLayout, 8> {

    using Base = DefaultLayoutImpl<CacheEfficientLayout, 8>;
    Index numElements;

    CacheEfficientLayout() noexcept : numElements(0) {}

    explicit CacheEfficientLayout(Index numElems) noexcept {
        ptr = makeUniqueForOverwrite<Elem>(allocatedSizeInElems(numElems));
        numElements = numElems;
        // TODO: ensure that allocated memory is cache aligned
    }

    constexpr static Index superblockSize() noexcept { return 4; }

    constexpr static Index blockSize() noexcept { return 1; }

    constexpr static Index bytesPerBlockCount() noexcept { return 1; }

    static Index allocatedSizeInElems(Index numElements) noexcept {
        return roundUpDiv(numElements, 4) * ELEMS_PER_CACHELINE;
    }

    [[nodiscard]] Index allocatedSizeInElems() const noexcept { return allocatedSizeInElems(numElems()); }

    [[nodiscard]] Index numElems() const noexcept { return numElements; }

    static Index getElemIdx(Index i) noexcept {
        assert(i >= 0);
        return (i / superblockSize()) * (ELEMS_PER_CACHELINE)
               + i % superblockSize(); // TODO: Make sure this gets optimized to bitshifts (possibly negative i may interfere)
    }

    [[nodiscard]] Elem& getElem(Index i) noexcept {
        assert(i >= 0 && i < numElements);
        return ptr[getElemIdx(i)];
    }
    [[nodiscard]] const Elem& getElem(Index i) const noexcept {
        assert(i >= 0 && i < numElements);
        return ptr[getElemIdx(i)];
    }

    void setElem(Index i, Elem value) noexcept {
        assert(i >= 0 && i < numElements);
        ptr[getElemIdx(i)] = value;
    }

    Elem& getElemRef(Index i) noexcept {
        assert(i >= 0 && i < numElements);
        return ptr[getElemIdx(i)];
    }

    [[nodiscard]] bool getBit(Index i) const noexcept {
        Elem v = getElem(i / 64);
        return BitwiseAccess<1>::getBits(&v, i % 64);
    }

    void setBit(Index i, bool value) noexcept { BitwiseAccess<1>::setBits(&getElemRef(i / 64), i % 64, value); }

    static Index getSuperblockCountIdx(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0);
        return superblockIdx * ELEMS_PER_CACHELINE + superblockSize();
    }

    [[nodiscard]] Elem getSuperblockCount(Index superblockIdx) const noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks());
        return ptr[getSuperblockCountIdx(superblockIdx)];
    }

    void setSuperblockCount(Index superblockIdx, Elem count) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks());
        ptr[getSuperblockCountIdx(superblockIdx)] = count;
    }

    [[nodiscard]] Index getBlockCount(Index blockIdx) const noexcept {
        Index superblockIdx = blockIdx / numBlocksInSuperblock();
        blockIdx %= numBlocksInSuperblock();
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks());
        return (Index)BitwiseAccess<8>::getBits(ptr.get(), getSuperblockCount(superblockIdx) + 1, blockIdx);
    }

    void setBlockCount(Index superblockIdx, Index blockIdx, Index value) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks() && blockIdx >= 0 && blockIdx < numBlocks());
        assert(value >= 0);
        Index i = getSuperblockCountIdx(superblockIdx) + 1;
        BitwiseAccess<8>::setBits(ptr.get(), i, blockIdx, value);
    }
};


/// TODO: Doxygen
// TODO: Test with other block sizes, fix bugs
template<Index SuperblockSize = (1 << 16) / 64, Index BlockSize = 8>
struct SimpleLayout : DefaultLayoutImpl<SimpleLayout<SuperblockSize, BlockSize>, BlockSize * 8> {

    static_assert(0 < BlockSize && 0 < SuperblockSize && SuperblockSize % BlockSize == 0);

private:
    ADS_CONSTEVAL static Index numBytesOfBlockCount() noexcept { return bytesNeededForIndexing(SuperblockSize * 64); }

public:
    using Base = DefaultLayoutImpl<SimpleLayout<SuperblockSize, BlockSize>, BlockSize * 8>;
    using BlockCount = IntType<numBytesOfBlockCount()>;
    using BlocksView = View<BlockCount>;
    Array<Elem> vec;
    View<Elem> superblocks;
    BlocksView blocks;

    using Base::numBlocks;
    using Base::numBlocksInSuperblock;
    using Base::numSuperblocks;

    SimpleLayout() noexcept = default;

    explicit SimpleLayout(Index numElements) : vec{makeUniqueForOverwrite<Elem>(allocatedSizeInElems(numElements))} {
        superblocks = View<Elem>(vec.ptr.get() + numElements, numSuperblocks(numElements));
        blocks = BlocksView(superblocks.ptr + numSuperblocks(), numBlocks(numElements));
    }

    [[nodiscard]] constexpr static Index superblockSize() noexcept { return SuperblockSize; }

    [[nodiscard]] constexpr static Index blockSize() noexcept { return BlockSize; }

    // TODO: Move to base class?
    [[nodiscard]] constexpr static Index bytesPerBlockCount() noexcept {
        constexpr Index res = numBytesOfBlockCount(numBlocksInSuperblock()); // force compile time evaluation in C++17 mode
        return res;
    }

    [[nodiscard]] Index allocatedSizeInElems() const noexcept { return allocatedSizeInElems(numElems()); }

    [[nodiscard]] static Index allocatedSizeInElems(Index numElements) noexcept {
        static_assert(sizeof(Elem) % sizeof(BlockCount) == 0);
        return numElements + blockSize() * numBlocks(numElements) * roundUpDiv(sizeof(Elem), sizeof(BlockCount))
               + superblockSize() * numSuperblocks(numElements);
    }

    [[nodiscard]] bool getBit(Index i) const noexcept { return bool(BitwiseAccess<1>::getBits(vec.ptr.get(), i)); }

    void setBit(Index i, bool value) noexcept { BitwiseAccess<1>::setBits(vec.ptr.get(), i, value); }

    [[nodiscard]] Elem& getElemRef(Index i) noexcept { return vec.bitAccess.getRef(vec.ptr.get(), i); }

    [[nodiscard]] Elem getElem(Index i) const noexcept { return vec[i]; }

    void setElem(Index i, Elem value) noexcept { vec.setBits(i, value); }

    [[nodiscard]] Elem getSuperblockCount(Index i) const noexcept { return superblocks[i]; }

    void setSuperblockCount(Index i, Elem value) noexcept { superblocks.setBits(i, value); }

    [[nodiscard]] Index getBlockCount(Index blockIdx) const noexcept { return blocks[blockIdx]; }

    void setBlockCount(Index superblockIdx, Index blockIdx, Index value) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks() && blockIdx >= 0 && blockIdx < numBlocksInSuperblock());
        assert(value >= 0);
        Index idx = superblockIdx * numBlocksInSuperblock() + blockIdx;
        blocks.setBits(idx, value);
    }

    [[nodiscard]] Index numElems() const noexcept { return superblocks.ptr - vec.ptr.get(); }
};

#ifdef ADS_HAS_CPP20
static_assert(IsLayout<CacheEfficientLayout>);
static_assert(IsLayout<SimpleLayout<>>);
static_assert(IsLayout<SimpleLayout<16>>);
#endif

} // namespace ads

#endif // BITVECTOR_BITVEC_LAYOUTS_HPP
