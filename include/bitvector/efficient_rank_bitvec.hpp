#ifndef ADS_EFFICIENT_RANK_BITVEC_HPP
#define ADS_EFFICIENT_RANK_BITVEC_HPP

#include "../bit_access.hpp"
#include "../common.hpp"
#include "bitvec_base.hpp"


namespace ads {
/// \brief An optimized Bitvector that can answer rank queries in O(1) but needs O(log n) for select queries.
/// \tparam BlockSizeInLimbs The number of Limbs in a block. This is the most important hyperparameter as it matters most for additional space requirements
///  and for performance. On at least some systems, the number of popcount invocations, which is in [0, BlockSize), is the most significant
/// factor for rank() performance.
/// \tparam SuperblockSizeInLimbs The number of Limbs in a superblock. This determines the amount of bytes needed to store a single block count, so this
/// parameter shouldn't be chosen too large. Must be a multiple of BlockSize.
/// \tparam SuperblockRankT The type used to store superblock counts. Some small memory requirement reductions are possible
/// if the size of the bitvector in bits is known in advance to fit into a smaller type than Limb.
template<Index BlockSizeInLimbs = 8, Index SuperblockSizeInLimbs = (1 << 16) / 64, typename SuperblockRankT = Limb>
class [[nodiscard]] EfficientRankBitvec
    : public BitvecBase<EfficientRankBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockRankT>,
              IntType<bytesNeededForIndexing(SuperblockSizeInLimbs * 64)>, SuperblockRankT> {

    using Base = BitvecBase<EfficientRankBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockRankT>,
            IntType<bytesNeededForIndexing(SuperblockSizeInLimbs * 64)>, SuperblockRankT>;
    friend Base;

    static_assert(0 < BlockSizeInLimbs && 0 < SuperblockSizeInLimbs && SuperblockSizeInLimbs % BlockSizeInLimbs == 0);


public:
    [[nodiscard]] constexpr static Index numLimbsInSuperblock() noexcept { return SuperblockSizeInLimbs; }

    [[nodiscard]] constexpr static Index numLimbsInBlock() noexcept { return BlockSizeInLimbs; }

    [[nodiscard]] ADS_CONSTEVAL static Index bytesPerBlockCount() noexcept {
        return bytesNeededForIndexing(Base::superblockSize());
    }

    using typename Base::BlockRank;
    static_assert(std::is_same_v<BlockRank, IntType<bytesPerBlockCount()>>);
    using BlockRanks = View<BlockRank>;
    using typename Base::SuperblockRank;
    static_assert(std::is_same_v<SuperblockRank, SuperblockRankT>);
    using SuperblockRanks = View<SuperblockRank>;
    Index length = 0;
    View<Limb> vec = View<Limb>();
    SuperblockRanks superblocks = SuperblockRanks();
    BlockRanks blocks = BlockRanks();

private:
    constexpr EfficientRankBitvec(UninitializedTag, Index numBits, Limb* mem = nullptr)
        : Base(numBits, mem), length(numBits) {
        Index numLimbs = roundUpDiv(numBits, 64);
        vec = View<Limb>(this->allocation.memory(), numLimbs);
        superblocks = SuperblockRanks(this->allocation.memory() + numLimbs, this->numSuperblocksForBits(numBits) + 1);
        blocks = BlockRanks(superblocks.ptr + this->numSuperblocks() + 1, this->numBlocksForBits(numBits));
    }

public:
    EfficientRankBitvec() noexcept = default;

    explicit constexpr EfficientRankBitvec(Index numBits, Limb fill, Limb* mem = nullptr)
        : EfficientRankBitvec(UninitializedTag{}, numBits, mem) {
        this->fill(fill);
    }

    explicit ADS_CPP20_CONSTEXPR EfficientRankBitvec(Span<const Limb> limbs, Limb* mem = nullptr) noexcept
        : EfficientRankBitvec(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR EfficientRankBitvec(Span<const Limb> limbs, Index numBits, Limb* mem = nullptr) noexcept
        : EfficientRankBitvec(UninitializedTag{}, numBits, mem) {
        this->copyFrom(limbs);
    }


    explicit ADS_CPP20_CONSTEXPR EfficientRankBitvec(std::string_view str, Index base = 2, Limb* mem = nullptr) noexcept
        : EfficientRankBitvec(UninitializedTag{}, str.size() * intLog2(base), mem) {
        this->initFromStr(str, base);
    }

    [[nodiscard]] static constexpr Index allocatedSizeInLimbsForLimbs(Index numLimbs) noexcept {
        static_assert(sizeof(Limb) % sizeof(BlockRank) == 0);
        static_assert(sizeof(Limb) % sizeof(SuperblockRank) == 0);
        Index numBits = numLimbs * 64;
        return numLimbs + roundUpDiv(Base::numBlocksForBits(numBits) * sizeof(BlockRank), sizeof(Limb))
               + roundUpDiv((1 + Base::numSuperblocksForBits(numBits)) * sizeof(SuperblockRank), sizeof(Limb));
    }

    [[nodiscard]] constexpr Index numLimbs() const noexcept {
        ADS_ASSUME(vec.numT == roundUpDiv(length, 64));
        return vec.numT;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numSuperblocks() const noexcept {
        ADS_ASSUME(superblocks.numT - 1 == roundUpDiv(length, this->superblockSize()));
        return superblocks.numT - 1; // the number of superblock rank counts is one larger
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBlocks() const noexcept {
        ADS_ASSUME(blocks.numT == roundUpDiv(length, this->blockSize()));
        return blocks.numT;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return length; }

private:
    [[nodiscard]] ADS_CPP20_CONSTEXPR const Limb* getLimbArray() const noexcept { return vec.ptr; }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb* getLimbArray() noexcept { return vec.ptr; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const View<BlockRank>& getBlockRankArray() const noexcept { return blocks; }
    [[nodiscard]] ADS_CPP20_CONSTEXPR View<BlockRank>& getBlockRankArray() noexcept { return blocks; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const View<SuperblockRank>& getSuperblockRankArray() const noexcept {
        return superblocks;
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR View<SuperblockRank>& getSuperblockRankArray() noexcept { return superblocks; }
};

template<typename = void>
using EfficientSelectBitvec = EfficientRankBitvec<>; // TODO: Remove
} // namespace ads

#endif // ADS_EFFICIENT_RANK_BITVEC_HPP
