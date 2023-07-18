#ifndef ADS_CLASSICAL_RANK_BITVEC_HPP
#define ADS_CLASSICAL_RANK_BITVEC_HPP

#include "../bit_access.hpp"
#include "../common.hpp"
#include "base/superblock_bitvec.hpp"


namespace ads {
/// \brief An optimized Bitvector that can answer rank queries in O(1) but needs O(log n) for select queries.
/// \tparam BlockSizeInLimbs The number of Limbs in a block. This is the most important hyperparameter as it matters most for additional space requirements
///  and for performance. On at least some systems, the number of popcount invocations, which is in [0, BlockSize), is the most significant
/// factor for rank() performance.
/// \tparam SuperblockSizeInLimbs The number of Limbs in a superblock. This determines the amount of bytes needed to store a single block count, so this
/// parameter shouldn't be chosen too large. Must be a multiple of BlockSize.
/// \tparam SuperblockRankT The type used to store superblock counts. Some small memory requirement reductions are possible
/// if the size of the bitvector in bits is known in advance to fit into a smaller type than Limb.
template<typename Derived, Index BlockSizeInLimbs, Index SuperblockSizeInLimbs, typename SuperblockRankT = Limb, typename OverwriteBlockRankT = FalseT>
class [[nodiscard]] ClassicalRankBitvecImpl
    : public SuperblockBitvecBase<Derived, SuperblockSizeInLimbs, BlockSizeInLimbs, SuperblockRankT, U64_PER_CACHELINE, OverwriteBlockRankT> {
protected:
    using Base = SuperblockBitvecBase<Derived, SuperblockSizeInLimbs, BlockSizeInLimbs, SuperblockRankT, U64_PER_CACHELINE, OverwriteBlockRankT>;
    friend Derived;
    friend Base;
    friend typename Base::Base;
    friend typename Base::Base::Base; // :D

    static_assert(0 < BlockSizeInLimbs && 0 < SuperblockSizeInLimbs && SuperblockSizeInLimbs % BlockSizeInLimbs == 0);


public:
    [[nodiscard]] constexpr static Index numLimbsInSuperblock() noexcept { return SuperblockSizeInLimbs; }

    [[nodiscard]] constexpr static Index numLimbsInBlock() noexcept { return BlockSizeInLimbs; }

    [[nodiscard]] ADS_CONSTEVAL static Index bytesPerBlockCount() noexcept {
        return bytesNeededForIndexing(SuperblockSizeInLimbs * 64);
    }

    using typename Base::BlockRank;
    static_assert(std::is_same_v<BlockRank, IntType<bytesPerBlockCount()>> || OverwriteBlockRankT{});
    using BlockRanks = Array<BlockRank>;
    using typename Base::SuperblockRank;
    static_assert(std::is_same_v<SuperblockRank, SuperblockRankT>);
    using SuperblockRanks = Array<SuperblockRank>;
    Array<CacheLine> vec = Array<CacheLine>();
    SuperblockRanks superblocks = SuperblockRanks();
    BlockRanks blocks = BlockRanks();

protected:
    [[nodiscard]] ADS_CPP20_CONSTEXPR const BlockRank* endOfMemory() const noexcept { return blocks.end(); }
    [[nodiscard]] ADS_CPP20_CONSTEXPR BlockRank* endOfMemory() noexcept { return blocks.end(); }


    // ctor is private but the direct base class is a friend, so an indirectly derived class can't call this.
    constexpr ClassicalRankBitvecImpl(UninitializedTag, Index numBits, CacheLine* mem = nullptr) noexcept
        : Base(numBits, mem) {
        // don't store incomplete blocks or cache lines
        numBits = roundUpTo(numBits, CACHELINE_SIZE_BYTES * 8);
        mem = this->allocation.memory();
        ADS_ASSUME_ALIGNED(mem, CACHELINE_SIZE_BYTES);
        Index numCacheLines = numBits / (CACHELINE_SIZE_BYTES * 8);
        ADS_ASSUME(numCacheLines * CACHELINE_SIZE_BYTES * 8 == numBits);
        vec = Array<CacheLine>(mem, numCacheLines);
        mem += numCacheLines;
        superblocks = SuperblockRanks(mem, this->numSuperblocksForBits(numBits) + 1);
        blocks = BlockRanks(superblocks.end(), numBits / Base::blockSize());
        ADS_ASSUME(this->allocation.sizeInBytes() >= this->derived().allocatedSizeInBytesForBits(numBits));
        ADS_IF_CONSTEVAL {}
        else {
            ADS_ASSUME((const char*)endOfMemory() - (const char*)this->allocation.memory() <= allocatedSizeInBytesForBits(numBits));
            ADS_ASSUME(allocatedSizeInBytesForBits(numBits)
                               - ((const char*)endOfMemory() - (const char*)this->allocation.memory())
                       <= CACHELINE_SIZE_BYTES);
        }
    }

public:
    ClassicalRankBitvecImpl() noexcept = default;

    explicit constexpr ClassicalRankBitvecImpl(Index numBits, Limb fill, CacheLine* mem = nullptr)
        : ClassicalRankBitvecImpl(UninitializedTag{}, numBits, mem) {
        this->fill(fill);
    }

    explicit ADS_CPP20_CONSTEXPR ClassicalRankBitvecImpl(Span<const Limb> limbs, CacheLine* mem = nullptr) noexcept
        : ClassicalRankBitvecImpl(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR ClassicalRankBitvecImpl(Span<const Limb> limbs, Index numBits, CacheLine* mem = nullptr) noexcept
        : ClassicalRankBitvecImpl(UninitializedTag{}, numBits, mem) {
        this->copyFrom(limbs);
    }


    explicit ADS_CPP20_CONSTEXPR ClassicalRankBitvecImpl(std::string_view str, Index base = 2, CacheLine* mem = nullptr) noexcept
        : ClassicalRankBitvecImpl(UninitializedTag{}, str.size() * intLog2(base), mem) {
        this->initFromStr(str, base);
    }

    [[nodiscard]] static constexpr Index allocatedSizeInBytesForBits(Index numBits) noexcept {
        numBits = roundUpTo(numBits, CACHELINE_SIZE_BYTES * 8);
        Index numCacheLines = roundUpDiv(numBits, 8 * CACHELINE_SIZE_BYTES);
        static_assert(sizeof(Limb) % sizeof(BlockRank) == 0);
        static_assert(sizeof(Limb) % sizeof(SuperblockRank) == 0);
        Index bytesForMetadata = numBits / Base::blockSize() * sizeof(BlockRank)
                                 + (1 + Base::numSuperblocksForBits(numBits)) * sizeof(SuperblockRank);
        return CACHELINE_SIZE_BYTES * numCacheLines + roundUpTo(bytesForMetadata, CACHELINE_SIZE_BYTES);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numSuperblocks() const noexcept {
        ADS_ASSUME(superblocks.numT - 1 == roundUpDiv(this->size(), this->superblockSize()));
        return superblocks.numT - 1; // the number of superblock rank counts is one larger
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numAccessibleBlocks() const noexcept {
        // the last block is always completed by appending zeros, so all blocks have the same size
        ADS_ASSUME(blocks.numT >= roundUpDiv(this->size(), this->blockSize()));
        ADS_ASSUME(blocks.numT * numLimbsInBlock() == this->numAccessibleLimbs());
        return blocks.numT;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<const CacheLine> getCacheLineArray() const noexcept {
        return {vec.ptr, vec.numT};
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<CacheLine> getCacheLineArray() noexcept { return {vec.ptr, vec.numT}; }


    [[nodiscard]] ADS_CPP20_CONSTEXPR BlockRank getBlockRank(Index i) const noexcept { return blocks[i]; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR SuperblockRank getSuperblockRank(Index i) const noexcept {
        return superblocks[i];
    }

    ADS_CPP20_CONSTEXPR void setBlockRank(Index i, BlockRank newVal) noexcept { blocks.setBits(i, newVal); }

    ADS_CPP20_CONSTEXPR void setBlockRank_(Index superblockIdx, Index blockInSuperblock, BlockRank newVal) noexcept {
        setBlockRank(superblockIdx * this->numBlocksInSuperblock() + blockInSuperblock, newVal);
    }

    ADS_CPP20_CONSTEXPR void setSuperblockRank(Index i, SuperblockRank newVal) noexcept {
        superblocks.setBits(i, newVal);
    }
};


template<Index BlockSizeInLimbs = 8, Index SuperblockSizeInLimbs = (1 << 16) / 64, typename SuperblockRankT = Limb>
class [[nodiscard]] ClassicalRankBitvec
    : public ClassicalRankBitvecImpl<ClassicalRankBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockRankT>, BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockRankT> {
public:
    using Base = ClassicalRankBitvecImpl<ClassicalRankBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockRankT>,
            BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockRankT>;
    using Base::Base;
    friend Base;

    constexpr ClassicalRankBitvec(UninitializedTag, Index numBits, CacheLine* mem = nullptr) noexcept
        : Base(UninitializedTag{}, numBits, mem) {}
};

#ifdef ADS_HAS_CPP20
static_assert(IsNormalBitvec<ClassicalRankBitvec<>>);
#endif

} // namespace ads

#endif // ADS_CLASSICAL_RANK_BITVEC_HPP
