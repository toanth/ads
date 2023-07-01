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
    { T::numElemsInSuperBlock() } -> std::convertible_to<Index>;
    { T::blockSize() } -> std::convertible_to<Index>;
    { T::numElemsInBlock() } -> std::convertible_to<Index>;
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


    using BlockCount = Elem;
    using SuperBlockCount = Elem;

    Allocation<> alloc; // must be the first member to prevent UB in the dtor
    View<Elem> vec = View<Elem>();
    Index numElements;

    constexpr CacheEfficientLayout() noexcept : numElements(0) {}

    explicit ADS_CPP20_CONSTEXPR CacheEfficientLayout(Index numElems, Elem* mem = nullptr) noexcept
        : alloc(allocatedSizeInElems(numElems), mem), vec(alloc.memory(), alloc.size()), numElements(numElems) {
        // TODO: ensure that allocated memory is cache aligned
    }


    constexpr static Index numElemsInSuperBlock() noexcept { return 4; }

    constexpr static Index superBlockSize() noexcept { return numElemsInSuperBlock() * 64; }

    constexpr static Index numElemsInBlock() noexcept { return 1; }

    constexpr static Index blockSize() noexcept { return numElemsInBlock() * 64; }

    constexpr static Index bytesPerBlockCount() noexcept { return 1; }


    [[nodiscard]] constexpr static Index numBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, numElemsInBlock());
    }

    [[nodiscard]] constexpr Index numBlocks() const noexcept { return numBlocks(numElems()); }

    [[nodiscard]] constexpr static Index numSuperBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, numElemsInSuperBlock());
    }

    [[nodiscard]] constexpr Index numSuperBlocks() const noexcept { return numSuperBlocks(numElems()); }

    [[nodiscard]] ADS_CONSTEVAL static Index numBlocksInSuperBlock() noexcept {
        assert(superBlockSize() % blockSize() == 0);
        return superBlockSize() / blockSize();
    }

    static constexpr Index allocatedSizeInElems(Index numElements) noexcept {
        return (1 + roundUpDiv(numElements, numElemsInSuperBlock())) * ELEMS_PER_CACHELINE;
    }

    [[nodiscard]] static constexpr Index allocatedSizeInElemsForBits(Index numBits) noexcept {
        return allocatedSizeInElems(roundUpDiv(numBits, 64));
    }

    [[nodiscard]] constexpr Index allocatedSizeInElems() const noexcept { return allocatedSizeInElems(numElems()); }

    [[nodiscard]] constexpr Index numElems() const noexcept { return numElements; }

    static constexpr Index getElemIdx(Index i) noexcept {
        ADS_ASSUME(i >= 0);
        return (i / numElemsInSuperBlock()) * (ELEMS_PER_CACHELINE) + i % numElemsInSuperBlock();
    }

    [[nodiscard]] constexpr Elem& getElem(Index i) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < numElements);
        return vec.ptr[i];
    }
    [[nodiscard]] constexpr Elem getElem(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < numElements);
        return vec[getElemIdx(i)];
    }

    constexpr void setElem(Index i, Elem value) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < numElements);
        vec.ptr[getElemIdx(i)] = value;
    }

    [[nodiscard]] constexpr Elem& getElemRef(Index i) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < numElements);
        return vec.ptr[getElemIdx(i)];
    }

    [[nodiscard]] constexpr bool getBit(Index i) const noexcept {
        Elem v = getElem(i / 64);
        return BitwiseAccess<1>::getBits(&v, i % 64);
    }

    constexpr void setBit(Index i, bool value) noexcept {
        BitwiseAccess<1>::setBits(&getElemRef(i / 64), i % 64, value);
    }

    [[nodiscard]] static constexpr Index getSuperBlockCountIdx(Index superblockIdx) noexcept {
        ADS_ASSUME(superblockIdx >= 0);
        return superblockIdx * ELEMS_PER_CACHELINE + numElemsInSuperBlock();
    }

    [[nodiscard]] constexpr Elem getSuperBlockCount(Index superblockIdx) const noexcept {
        return vec[getSuperBlockCountIdx(superblockIdx)];
    }

    constexpr void setSuperBlockCount(Index superblockIdx, Elem count) noexcept {
        vec.ptr[getSuperBlockCountIdx(superblockIdx)] = count;
    }

    [[nodiscard]] constexpr Index getBlockCount(Index blockIdx) const noexcept {
        Index superblockIdx = blockIdx / numBlocksInSuperBlock();
        blockIdx %= numBlocksInSuperBlock();
        ADS_ASSUME(superblockIdx >= 0);
        ADS_ASSUME(superblockIdx < numSuperBlocks());
        return (Index)BitwiseAccess<8>::getBits(vec.ptr, getSuperBlockCountIdx(superblockIdx) + 1, blockIdx);
    }

    constexpr void setBlockCount(Index superblockIdx, Index blockIdx, Index value) noexcept {
        ADS_ASSUME(superblockIdx >= 0 && superblockIdx < numSuperBlocks() && blockIdx >= 0 && blockIdx < numBlocks());
        ADS_ASSUME(value >= 0);
        Index i = getSuperBlockCountIdx(superblockIdx) + 1;
        BitwiseAccess<8>::setBits(vec.ptr, i, blockIdx, value);
    }
};


/// \brief The default layout for the Bitvector.
/// \tparam BlockSizeInElems The number of Elems in a block. This is the most important hyperparameter as it matters most for additional space requirements
///  and for performance. On at least some systems, the number of popcount invocations, which is in [0, BlockSize), is the most significant
/// factor for rank() performance.
/// \tparam SuperBlockSizeInElems The number of Elems in a superblock. This determines the amount of bytes needed to store a single block count, so this
/// parameter shouldn't be chosen too large. Must be a multiple of BlockSize.
/// \tparam SuperBlockCountT The type used to store superblock counts. Some small memory requirement reductions are possible
/// if the size of the bitvector in bits is known in advance to fit into a smaller type than Elem.
template<Index BlockSizeInElems = 8, Index SuperBlockSizeInElems = (1 << 16) / 64, typename SuperBlockCountT = Elem>
struct SimpleLayout {

    static_assert(0 < BlockSizeInElems && 0 < SuperBlockSizeInElems && SuperBlockSizeInElems % BlockSizeInElems == 0);


public:
    [[nodiscard]] constexpr static Index numElemsInSuperBlock() noexcept { return SuperBlockSizeInElems; }
    [[nodiscard]] constexpr static Index superBlockSize() noexcept { return SuperBlockSizeInElems * 64; }

    [[nodiscard]] constexpr static Index numElemsInBlock() noexcept { return BlockSizeInElems; }
    [[nodiscard]] constexpr static Index blockSize() noexcept { return BlockSizeInElems * 64; }

    [[nodiscard]] ADS_CONSTEVAL static Index bytesPerBlockCount() noexcept {
        return bytesNeededForIndexing(superBlockSize());
    }

    using BlockCount = IntType<bytesPerBlockCount()>;
    using BlockCounts = View<BlockCount>;
    using SuperBlockCount = SuperBlockCountT;
    using SuperBlockCounts = View<SuperBlockCount>;
    Allocation<> alloc; // must be the first data member to prevent UB in the dtor
    View<Elem> vec;
    SuperBlockCounts superBlocks;
    BlockCounts blocks;

    SimpleLayout() noexcept = default;

    explicit constexpr SimpleLayout(Index numElements, Elem* mem = nullptr)
        : alloc(allocatedSizeInElems(numElements), mem), vec(alloc.memory(), numElements) {
        superBlocks = SuperBlockCounts(alloc.memory() + numElements, numSuperBlocks(numElements) + 1);
        blocks = BlockCounts(superBlocks.ptr + numSuperBlocks() + 1, numBlocks(numElements));
    }

    [[nodiscard]] constexpr static Index numBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, numElemsInBlock());
    }

    [[nodiscard]] constexpr Index numBlocks() const noexcept { return numBlocks(numElems()); }

    [[nodiscard]] constexpr static Index numSuperBlocks(Index numElems) noexcept {
        return roundUpDiv(numElems, numElemsInSuperBlock());
    }

    [[nodiscard]] constexpr Index numSuperBlocks() const noexcept { return numSuperBlocks(numElems()); }

    [[nodiscard]] ADS_CONSTEVAL static Index numBlocksInSuperBlock() noexcept {
        assert(superBlockSize() % blockSize() == 0);
        return superBlockSize() / blockSize();
    }

    [[nodiscard]] constexpr Index allocatedSizeInElems() const noexcept { return allocatedSizeInElems(numElems()); }

    [[nodiscard]] static constexpr Index allocatedSizeInElems(Index numElements) noexcept {
        static_assert(sizeof(Elem) % sizeof(BlockCount) == 0);
        static_assert(sizeof(Elem) % sizeof(SuperBlockCount) == 0);
        return numElements + roundUpDiv(numBlocks(numElements) * sizeof(BlockCount), sizeof(Elem))
               + roundUpDiv((1 + numSuperBlocks(numElements)) * sizeof(SuperBlockCount), sizeof(Elem));
    }

    [[nodiscard]] static constexpr Index allocatedSizeInElemsForBits(Index numBits) noexcept {
        return allocatedSizeInElems(roundUpDiv(numBits, 64));
    }

    [[nodiscard]] constexpr bool getBit(Index i) const noexcept { return bool(BitwiseAccess<1>::getBits(vec.ptr, i)); }

    constexpr void setBit(Index i, bool value) noexcept { BitwiseAccess<1>::setBits(vec.ptr, i, value); }

    [[nodiscard]] constexpr Elem& getElemRef(Index i) noexcept { return vec.bitAccess.getRef(vec.ptr, i); }

    [[nodiscard]] constexpr Elem getElem(Index i) const noexcept { return vec[i]; }

    constexpr void setElem(Index i, Elem value) noexcept { vec.setBits(i, value); }

    [[nodiscard]] constexpr Elem getSuperBlockCount(Index i) const noexcept { return superBlocks[i]; }

    constexpr void setSuperBlockCount(Index i, Elem value) noexcept { superBlocks.setBits(i, value); }

    [[nodiscard]] constexpr Index getBlockCount(Index blockIdx) const noexcept { return blocks[blockIdx]; }

    constexpr void setBlockCount(Index superBlockIdx, Index blockIdx, Index value) noexcept {
        assert(superBlockIdx >= 0 && superBlockIdx < numSuperBlocks() && blockIdx >= 0 && blockIdx < numBlocksInSuperBlock());
        assert(value >= 0);
        Index idx = superBlockIdx * numBlocksInSuperBlock() + blockIdx;
        blocks.setBits(idx, value);
    }

    [[nodiscard]] constexpr Index numElems() const noexcept { return vec.numT; }
};

#ifdef ADS_HAS_CPP20
static_assert(IsLayout<CacheEfficientLayout>);
static_assert(IsLayout<SimpleLayout<1>>);
static_assert(IsLayout<SimpleLayout<16>>);
#endif

} // namespace ads

#endif // BITVECTOR_BITVEC_LAYOUTS_HPP
