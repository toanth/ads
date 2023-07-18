#ifndef ADS_CLASSICAL_BITVEC_HPP
#define ADS_CLASSICAL_BITVEC_HPP

#include "classical_rank_bitvec.hpp"

namespace ads {

template<typename Derived, Index BlockSizeInLimbs, Index SuperblockSizeInLimbs,
        Index SuperblockSelectStepSize = SuperblockSizeInLimbs * 64, Index InSuperblockSelectStepSize = BlockSizeInLimbs * 64,
        typename SuperblockSelectT = U32, typename SuperblockRankT = Limb, typename OverwriteRankT = FalseT>
class ClassicalBitvecImpl
    : public ClassicalRankBitvecImpl<Derived, BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockRankT, OverwriteRankT> {
protected:
    using Base = ClassicalRankBitvecImpl<Derived, BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockRankT, OverwriteRankT>;
    using Base::Base;
    friend Base;

    using typename Base::BlockRank;
    using InSuperblockSelectT = IntType<bytesNeededForIndexing(SuperblockSizeInLimbs / BlockSizeInLimbs)>;

    static_assert((SuperblockSizeInLimbs * 64) % InSuperblockSelectStepSize == 0);
    constexpr static Index numSelectBlocksInSuperblock = SuperblockSizeInLimbs * 64 / InSuperblockSelectStepSize;

    /// \brief The ith value stores the index of the superblocks where the SuperblockSelectStepSize * ith one is and
    /// the size - 1 - ith index stores the index of the SuperblockSelectStepSize * ith zero.
    Array<SuperblockSelectT> superblockSelect;
    /// \brief For each superblock i, store select one starting at i * numSelectBlocksInSuperblock and select zero
    /// starting at (i + 1) * numSelectBlockInSuperblock - 1, growing downwards
    Array<InSuperblockSelectT> selectInSuperblock;


public:
    constexpr ClassicalBitvecImpl(UninitializedTag, Index numBits, CacheLine* mem = nullptr) noexcept
        : Base(UninitializedTag{}, numBits, mem) {
        numBits = roundUpTo(numBits, CACHELINE_SIZE_BYTES * 8);
        Index numSuperblockSelectEntries = roundUpDiv(numBits, SuperblockSelectStepSize) + 2;
        ADS_ASSUME(Base::allocatedSizeInBytesForBits(numBits) % CACHELINE_SIZE_BYTES == 0);
        mem = this->allocation.memory() + Base::allocatedSizeInBytesForBits(numBits) / CACHELINE_SIZE_BYTES;
        superblockSelect = Array<SuperblockSelectT>(mem, numSuperblockSelectEntries);
        Index numSuperblocks = roundUpDiv(numBits, Base::superblockSize());
        selectInSuperblock
                = Array<InSuperblockSelectT>(superblockSelect.end(), numSuperblocks * (numSelectBlocksInSuperblock + 2));
        ADS_ASSUME(this->allocation.isEnd(endOfMemory()));
    }


    explicit constexpr ClassicalBitvecImpl(Index numBits, Limb fill, CacheLine* mem = nullptr)
        : ClassicalBitvecImpl(UninitializedTag{}, numBits, mem) {
        this->fill(fill);
    }

    explicit ADS_CPP20_CONSTEXPR ClassicalBitvecImpl(Span<const Limb> limbs, CacheLine* mem = nullptr) noexcept
        : ClassicalBitvecImpl(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR ClassicalBitvecImpl(Span<const Limb> limbs, Index numBits, CacheLine* mem = nullptr) noexcept
        : ClassicalBitvecImpl(UninitializedTag{}, numBits, mem) {
        this->copyFrom(limbs);
    }


    explicit ADS_CPP20_CONSTEXPR ClassicalBitvecImpl(std::string_view str, Index base = 2, CacheLine* mem = nullptr) noexcept
        : ClassicalBitvecImpl(UninitializedTag{}, str.size() * intLog2(base), mem) {
        this->initFromStr(str, base);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const InSuperblockSelectT* endOfMemory() const noexcept {
        return selectInSuperblock.end();
    }

    [[nodiscard]] static constexpr Index allocatedSizeInBytesForBits(Index numBits) noexcept {
        numBits = roundUpTo(numBits, CACHELINE_SIZE_BYTES * 8);
        Index baseSize = Base::allocatedSizeInBytesForBits(numBits);
        ADS_ASSUME(baseSize % CACHELINE_SIZE_BYTES == 0);
        Index numSuperblocks = roundUpDiv(numBits, Base::superblockSize());
        Index numSuperblockSelectEntries = roundUpDiv(numBits, SuperblockSelectStepSize) + 2;
        Index extra = sizeof(SuperblockSelectT) * numSuperblockSelectEntries
                      + sizeof(InSuperblockSelectT) * numSuperblocks * (numSelectBlocksInSuperblock + 2);
        return baseSize + roundUpTo(extra, CACHELINE_SIZE_BYTES);
    }


    ADS_CPP20_CONSTEXPR void buildSelectMetadata() noexcept {
        Index selectOneIdx = 0;
        Index selectZeroIdx = 0;
        for (Index sbIdx = 1; sbIdx < this->numSuperblocks() + 1; ++sbIdx) {
            Index rank = this->superblocks[sbIdx];
            while (rank > selectOneIdx * SuperblockSelectStepSize) {
                superblockSelect.setBits(selectOneIdx++, sbIdx - 1);
            }
            Index rankZero = sbIdx == this->numSuperblocks() ? this->size() - rank : sbIdx * this->superblockSize() - rank;
            while (rankZero > selectZeroIdx * SuperblockSelectStepSize) {
                superblockSelect.setBits(superblockSelect.numT - ++selectZeroIdx, sbIdx - 1);
            }
        }
        ADS_ASSUME(selectOneIdx <= superblockSelect.numT - 1 - selectZeroIdx);
        superblockSelect.setBits(selectOneIdx, this->numSuperblocks() - 1);
        superblockSelect.setBits(superblockSelect.numT - 1 - selectZeroIdx, this->numSuperblocks() - 1);

        for (Index sbIdx = 0; sbIdx < this->numSuperblocks(); ++sbIdx) {
            Index selectOneOffset = sbIdx * (numSelectBlocksInSuperblock + 2);
            Index selectZeroOffset = (sbIdx + 1) * (numSelectBlocksInSuperblock + 2) - 1;
            selectOneIdx = 0;
            selectZeroIdx = 0;
            Index numBlocks = (sbIdx + 1 == this->numSuperblocks() ? this->numBlocks() - sbIdx * this->numBlocksInSuperblock() :
                                                                     this->numBlocksInSuperblock());
            //            Index lastRank = 0;
            for (Index blockIdx = 0; blockIdx <= numBlocks; ++blockIdx) {
                Index rank = blockIdx == numBlocks ?
                                     this->superblocks[sbIdx + 1] :
                                     this->superblocks[sbIdx] + this->blocks[sbIdx * this->numBlocksInSuperblock() + blockIdx];
                //                if constexpr (OverwriteRankT{}) {
                //                    if (rank < lastRank || (rank == lastRank && this->getLimb(blockIdx * this->numLimbsInBlock()) != 0)) {
                //                    }
                //                }
                //                Index lastRank = rank;
                rank -= this->superblocks[sbIdx];
                Index rankZero = blockIdx * this->blockSize() - rank;
                while (rank > selectOneIdx * InSuperblockSelectStepSize) {
                    selectInSuperblock.setBits(selectOneOffset + selectOneIdx++, blockIdx - 1);
                }
                while (rankZero > selectZeroIdx * InSuperblockSelectStepSize) {
                    selectInSuperblock.setBits(selectZeroOffset - selectZeroIdx++, blockIdx - 1);
                }
                ADS_ASSUME(selectZeroOffset - selectZeroIdx >= selectOneOffset + selectOneIdx);
                selectInSuperblock.setBits(selectOneOffset + selectOneIdx, numBlocks - 1);
                selectInSuperblock.setBits(selectZeroOffset - selectZeroIdx, numBlocks - 1);
            }
        }
    }

    template<bool IsOne>
    [[nodiscard]] [[gnu::hot]] ADS_CPP20_CONSTEXPR Index selectSuperblockIdx(Index& bitRank) const noexcept {
        ADS_ASSUME(bitRank >= 0);
        ADS_ASSUME(bitRank < this->template numBitsEqualTo<IsOne>());
        Index inSuperblockSelectArr = bitRank / SuperblockSelectStepSize;
        auto getSuperblock = [this](Index superblockSelectIdx) {
            if constexpr (IsOne) {
                return superblockSelect[superblockSelectIdx];
            } else {
                return superblockSelect[superblockSelect.numT - 1 - superblockSelectIdx];
            }
        };
        Index superblockLowerBound = getSuperblock(inSuperblockSelectArr);
        // + 1 because we're searching for the first superblock where the rank is greater than bitRank
        //        ++superblockLowerBound; // TODO: Enable to improve performance?
        Index superblockUpperBound = getSuperblock(inSuperblockSelectArr + 1) + 1; // add 1 because the range is semi=open
        ADS_ASSUME(superblockLowerBound >= 0);
        ADS_ASSUME(superblockLowerBound < this->numSuperblocks());
        ADS_ASSUME(superblockLowerBound < superblockUpperBound);
        return this->template selectSuperblockIdxInRange<IsOne>(bitRank, superblockLowerBound, superblockUpperBound);
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectBlockIdx(Index& bitRank, Index superblockIdx) const noexcept {
        ADS_ASSUME(bitRank >= 0);
        ADS_ASSUME(bitRank < this->superblockSize());
        ADS_ASSUME(superblockIdx >= 0);
        ADS_ASSUME(superblockIdx < this->numSuperblocks());
        auto selectEntry = [](Index superblockIdx, Index inSuperblock) {
            ADS_ASSUME(inSuperblock >= 0);
            ADS_ASSUME(inSuperblock <= numSelectBlocksInSuperblock);
            if constexpr (IsOne) {
                return superblockIdx * (numSelectBlocksInSuperblock + 2) + inSuperblock;
            } else {
                Index endOfBlock = (superblockIdx + 1) * (numSelectBlocksInSuperblock + 2) - 1;
                return endOfBlock - inSuperblock;
            }
        };
        Index inSuperblock = bitRank / InSuperblockSelectStepSize;
        Index lowerSelectIdx = selectEntry(superblockIdx, inSuperblock);
        Index upperSelectIdx = selectEntry(superblockIdx, inSuperblock + 1);
        ADS_ASSUME(lowerSelectIdx >= 0);
        Index blockLowerBound = selectInSuperblock[lowerSelectIdx];
        Index blockUpperBound = selectInSuperblock[upperSelectIdx] + 1;
        ADS_ASSUME(blockLowerBound >= 0);
        ADS_ASSUME(blockUpperBound > blockLowerBound);
        ADS_ASSUME(blockUpperBound <= this->numBlocksInSuperblock());
        blockLowerBound += superblockIdx * this->numBlocksInSuperblock();
        blockUpperBound += superblockIdx * this->numBlocksInSuperblock();
        if constexpr (OverwriteRankT{}) {
            bitRank %= Limb(1) << (sizeof(Unwrap<OverwriteRankT>) * 8);
        }
        return this->template selectBlockIdxInRange<IsOne>(bitRank, superblockIdx, blockLowerBound, blockUpperBound);
    }
};


template<Index BlockSizeInLimbs = 8, Index SuperblockSizeInLimbs = (1 << 16) / 64 / 4, Index SuperblockSelectStepSize = SuperblockSizeInLimbs * 64,
        Index InSuperblockSelectStepSize = BlockSizeInLimbs * 64, typename SuperblockSelectT = U32, typename SuperblockRankT = Limb>
class ClassicalBitvec
    : public ClassicalBitvecImpl<ClassicalBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectStepSize, InSuperblockSelectStepSize, SuperblockSelectT, SuperblockRankT>,
              BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectStepSize, InSuperblockSelectStepSize, SuperblockSelectT, SuperblockRankT> {
public:
    using Base = ClassicalBitvecImpl<ClassicalBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectStepSize, InSuperblockSelectStepSize, SuperblockSelectT, SuperblockRankT>,
            BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectStepSize, InSuperblockSelectStepSize, SuperblockSelectT, SuperblockRankT>;
    using Base::Base;
    friend Base;

    constexpr ClassicalBitvec(UninitializedTag, Index numBits, CacheLine* mem = nullptr) noexcept
        : Base(UninitializedTag{}, numBits, mem) {}
};

} // namespace ads

#endif // ADS_CLASSICAL_BITVEC_HPP
