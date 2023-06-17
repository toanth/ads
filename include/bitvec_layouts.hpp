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
    { T::superBlockSize() } -> std::convertible_to<Index>;
    { T::blockSize() } -> std::convertible_to<Index>;
    { T::numBlocksInSuperBlock() } -> std::convertible_to<Index>;
    { T::bytesPerBlockCount() } -> std::convertible_to<Index>;
    { ct.allocatedSizeInElems() } -> std::convertible_to<Index>;
    { ct.getElem(Index()) } -> std::convertible_to<const Elem&>;
    { ct.getSuperBlockCount(Index()) } -> std::convertible_to<const Elem&>;
    t.setSuperBlockCount(Index(), Elem());
    { ct.getBlockCount(Index()) } -> std::convertible_to<Index>;
    t.setBlockCount(Index(), Index(), Index());
    { ct.numElems() } -> std::convertible_to<Index>;
    { ct.numBlocks() } -> std::convertible_to<Index>;
    { ct.numSuperBlocks() } -> std::convertible_to<Index>;
    { T::allocatedSizeInElems(Index()) } -> std::convertible_to<Index>;
};
#define ADS_LAYOUT_CONCEPT IsLayout
#else
#define ADS_LAYOUT_CONCEPT class
#endif


// This layout wastes a lot of memory; it needs roughly double the total memory of the actual bit vector
// Still, for small bitvectors, this may be a good option
struct CacheEfficientLayout {

    Array<Elem> vec = Array<Elem>();
    Index numElements;

    CacheEfficientLayout() noexcept : numElements(0) {}

    explicit CacheEfficientLayout(Index numElems) noexcept {
        vec.ptr = makeUniqueForOverwrite<Elem>(allocatedSizeInElems(numElems));
        numElements = numElems;
        // TODO: ensure that allocated memory is cache aligned
    }

    constexpr static Index superBlockSize() noexcept { return 4; }

    constexpr static Index blockSize() noexcept { return 1; }

    constexpr static Index bytesPerBlockCount() noexcept { return 1; }


    [[nodiscard]] constexpr static Index numBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, blockSize());
    }

    [[nodiscard]] constexpr Index numBlocks() const noexcept { return roundUpDiv(numElems(), blockSize()); }

    [[nodiscard]] constexpr static Index numSuperBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, superBlockSize());
    }

    [[nodiscard]] constexpr Index numSuperBlocks() const noexcept { return numSuperBlocks(numElems()); }

    [[nodiscard]] ADS_CONSTEVAL static Index numBlocksInSuperBlock() noexcept {
        assert(superBlockSize() % blockSize() == 0);
        return superBlockSize() / blockSize();
    }

    static constexpr Index allocatedSizeInElems(Index numElements) noexcept {
        return roundUpDiv(numElements, 4) * ELEMS_PER_CACHELINE;
    }

    [[nodiscard]] constexpr Index allocatedSizeInElems() const noexcept { return allocatedSizeInElems(numElems()); }

    [[nodiscard]] constexpr Index numElems() const noexcept { return numElements; }

    static constexpr Index getElemIdx(Index i) noexcept {
        assert(i >= 0);
        return (i / superBlockSize()) * (ELEMS_PER_CACHELINE)
               + i % superBlockSize(); // TODO: Make sure this gets optimized to bitshifts (possibly negative i may interfere)
    }

    [[nodiscard]] Elem& getElem(Index i) noexcept {
        assert(i >= 0 && i < numElements);
        return vec.ptr[i];
    }
    [[nodiscard]] Elem getElem(Index i) const noexcept {
        assert(i >= 0 && i < numElements);
        return vec[getElemIdx(i)];
    }

    void setElem(Index i, Elem value) noexcept {
        assert(i >= 0 && i < numElements);
        vec.ptr[getElemIdx(i)] = value;
    }

    Elem& getElemRef(Index i) noexcept {
        assert(i >= 0 && i < numElements);
        return vec.ptr[getElemIdx(i)];
    }

    [[nodiscard]] bool getBit(Index i) const noexcept {
        Elem v = getElem(i / 64);
        return BitwiseAccess<1>::getBits(&v, i % 64);
    }

    void setBit(Index i, bool value) noexcept { BitwiseAccess<1>::setBits(&getElemRef(i / 64), i % 64, value); }

    static Index getSuperblockCountIdx(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0);
        return superblockIdx * ELEMS_PER_CACHELINE + superBlockSize();
    }

    [[nodiscard]] Elem getSuperBlockCount(Index superblockIdx) const noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperBlocks());
        return vec[getSuperblockCountIdx(superblockIdx)];
    }

    void setSuperBlockCount(Index superblockIdx, Elem count) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperBlocks());
        vec.ptr[getSuperblockCountIdx(superblockIdx)] = count;
    }

    [[nodiscard]] Index getBlockCount(Index blockIdx) const noexcept {
        Index superblockIdx = blockIdx / numBlocksInSuperBlock();
        blockIdx %= numBlocksInSuperBlock();
        assert(superblockIdx >= 0 && superblockIdx < numSuperBlocks());
        return (Index)BitwiseAccess<8>::getBits(vec.ptr.get(), getSuperBlockCount(superblockIdx) + 1, blockIdx);
    }

    void setBlockCount(Index superblockIdx, Index blockIdx, Index value) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperBlocks() && blockIdx >= 0 && blockIdx < numBlocks());
        assert(value >= 0);
        Index i = getSuperblockCountIdx(superblockIdx) + 1;
        BitwiseAccess<8>::setBits(vec.ptr.get(), i, blockIdx, value);
    }
};


/// TODO: Doxygen
// TODO: Test with other block sizes, fix bugs
template<Index SuperblockSize = (1 << 16) / 64, Index BlockSize = 1, typename SuperBlockCount = Elem>
struct SimpleLayout {

    static_assert(0 < BlockSize && 0 < SuperblockSize && SuperblockSize % BlockSize == 0);


public:
    [[nodiscard]] constexpr static Index superBlockSize() noexcept { return SuperblockSize; }

    [[nodiscard]] constexpr static Index blockSize() noexcept { return BlockSize; }

    [[nodiscard]] ADS_CONSTEVAL static Index bytesPerBlockCount() noexcept {
        return bytesNeededForIndexing(SuperblockSize * 64);
    }

    using BlockCount = IntType<bytesPerBlockCount()>;
    using BlockCounts = View<BlockCount>;
    using SuperBlockCounts = View<SuperBlockCount>;
    Array<Elem> vec;
    BlockCounts blocks;
    SuperBlockCounts superBlocks;

    SimpleLayout() noexcept = default;

    explicit SimpleLayout(Index numElements) : vec{makeUniqueForOverwrite<Elem>(allocatedSizeInElems(numElements))} {
        superBlocks = SuperBlockCounts(vec.ptr.get() + numElements, numSuperBlocks(numElements));
        blocks = BlockCounts(superBlocks.ptr + numSuperBlocks(), numBlocks(numElements));
    }

    [[nodiscard]] constexpr static Index numBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, blockSize());
    }

    [[nodiscard]] constexpr Index numBlocks() const noexcept { return numBlocks(numElems()); }

    [[nodiscard]] constexpr static Index numSuperBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, superBlockSize());
    }

    [[nodiscard]] constexpr Index numSuperBlocks() const noexcept { return numSuperBlocks(numElems()); }

    [[nodiscard]] ADS_CONSTEVAL static Index numBlocksInSuperBlock() noexcept {
        assert(superBlockSize() % blockSize() == 0);
        return superBlockSize() / blockSize();
    }

    [[nodiscard]] Index allocatedSizeInElems() const noexcept { return allocatedSizeInElems(numElems()); }

    [[nodiscard]] static Index allocatedSizeInElems(Index numElements) noexcept {
        static_assert(sizeof(Elem) % sizeof(BlockCount) == 0);
        static_assert(sizeof(Elem) % sizeof(SuperBlockCount) == 0);
        return numElements + roundUpDiv(numBlocks(numElements) * sizeof(BlockCount), sizeof(Elem))
               + roundUpDiv(numSuperBlocks(numElements) * sizeof(SuperBlockCount), sizeof(Elem));
    }

    [[nodiscard]] bool getBit(Index i) const noexcept { return bool(BitwiseAccess<1>::getBits(vec.ptr.get(), i)); }

    void setBit(Index i, bool value) noexcept { BitwiseAccess<1>::setBits(vec.ptr.get(), i, value); }

    [[nodiscard]] Elem& getElemRef(Index i) noexcept { return vec.bitAccess.getRef(vec.ptr.get(), i); }

    [[nodiscard]] Elem getElem(Index i) const noexcept { return vec[i]; }

    void setElem(Index i, Elem value) noexcept { vec.setBits(i, value); }

    [[nodiscard]] Elem getSuperBlockCount(Index i) const noexcept { return superBlocks[i]; }

    void setSuperBlockCount(Index i, Elem value) noexcept { superBlocks.setBits(i, value); }

    [[nodiscard]] Index getBlockCount(Index blockIdx) const noexcept { return blocks[blockIdx]; }

    void setBlockCount(Index superBlockIdx, Index blockIdx, Index value) noexcept {
        assert(superBlockIdx >= 0 && superBlockIdx < numSuperBlocks() && blockIdx >= 0 && blockIdx < numBlocksInSuperBlock());
        assert(value >= 0);
        Index idx = superBlockIdx * numBlocksInSuperBlock() + blockIdx;
        blocks.setBits(idx, value);
    }

    [[nodiscard]] Index numElems() const noexcept {
        assert(((const char*)(superBlocks.ptr) - (const char*)vec.ptr.get()) % sizeof(Elem) == 0);
        return reinterpret_cast<const Elem*>(superBlocks.ptr) - vec.ptr.get();
    }
};

#ifdef ADS_HAS_CPP20
static_assert(IsLayout<CacheEfficientLayout>);
static_assert(IsLayout<SimpleLayout<>>);
static_assert(IsLayout<SimpleLayout<16>>);
#endif

} // namespace ads

#endif // BITVECTOR_BITVEC_LAYOUTS_HPP
