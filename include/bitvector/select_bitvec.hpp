#ifndef ADS_SELECT_BITVEC_HPP
#define ADS_SELECT_BITVEC_HPP

#include "../common.hpp"
#include "classical_bitvec.hpp"

namespace ads {

// template<Index BlockSizeInLimbs = 4, Index SuperblockSizeInLimbs = (1 << 16) / 64 / 4, typename SuperblockSelectT = U32>
// class SelectBitvec
//     : public ClassicalBitvecImpl<SelectBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectT>, BlockSizeInLimbs,
//               SuperblockSizeInLimbs, SuperblockSizeInLimbs * 64, BlockSizeInLimbs * 64, SuperblockSelectT, Limb> {
//     //              SuperblockSizeInLimbs * 64, BlockSizeInLimbs * 64, SuperblockSelectT, Limb, OverwriteType<Byte>> {
//     using Base = ClassicalBitvecImpl<SelectBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectT>, BlockSizeInLimbs,
//             SuperblockSizeInLimbs, SuperblockSizeInLimbs * 64, BlockSizeInLimbs * 64, SuperblockSelectT, Limb>;
//     using Base::Base;
//     friend Base;
//     static_assert(BlockSizeInLimbs <= 4, "BlockSize not supported");
//     static_assert(BlockSizeInLimbs * 256 >= SuperblockSizeInLimbs, "SuperblockSize not supported");
//
// public:
//     constexpr SelectBitvec(UninitializedTag, Index numBits, CacheLine* mem = nullptr) noexcept
//         : Base(UninitializedTag{}, numBits, mem) {}
//
//     [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOneUnchecked(Index pos) const noexcept {
//         return this->template rankFallback<true>(pos); // not fast, but that's not a problem
//     }
//
//     template<bool IsOne>
//     [[nodiscard]] [[using gnu: hot]] ADS_CPP20_CONSTEXPR Index selectBlockIdx(Index& bitRank, Index superblockIdx) const noexcept {
//         ADS_ASSUME(bitRank >= 0);
//         ADS_ASSUME(bitRank < this->superblockSize());
//         ADS_ASSUME(superblockIdx >= 0);
//         ADS_ASSUME(superblockIdx < this->numSuperblocks());
//         auto selectEntry = [](Index superblockIdx, Index inSuperblock) {
//             ADS_ASSUME(inSuperblock >= 0);
//             ADS_ASSUME(inSuperblock <= SelectBitvec::numBlocksInSuperblock());
//             if constexpr (IsOne) {
//                 return superblockIdx * (SelectBitvec::numBlocksInSuperblock() + 2) + inSuperblock;
//             } else {
//                 Index endOfBlock = (superblockIdx + 1) * (SelectBitvec::numBlocksInSuperblock() + 2) - 1;
//                 return endOfBlock - inSuperblock;
//             }
//         };
//         Index inSuperblock = bitRank / this->blockSize();
//         Index lowerSelectIdx = selectEntry(superblockIdx, inSuperblock);
//         Index upperSelectIdx = selectEntry(superblockIdx, inSuperblock + 1);
//         ADS_ASSUME(lowerSelectIdx >= 0);
//         ADS_ASSUME(lowerSelectIdx <= upperSelectIdx);
//         Index blockLowerBound = this->selectInSuperblock[lowerSelectIdx];
//         Index blockUpperBound = this->selectInSuperblock[upperSelectIdx] + 1;
//         ADS_ASSUME(blockLowerBound >= 0);
//         ADS_ASSUME(blockUpperBound > blockLowerBound);
//         ADS_ASSUME(blockUpperBound <= this->numBlocksInSuperblock());
//         Index inSuperblockSoFar = lowerSelectIdx * this->blockSize();
//         blockLowerBound += superblockIdx * this->numBlocksInSuperblock();
//         blockUpperBound += superblockIdx * this->numBlocksInSuperblock();
//         return inSuperblockSoFar + this->template selectBlockIdxInRange<IsOne>(bitRank, superblockIdx, blockLowerBound, blockUpperBound);
//     }
// };


// template<Index BlockSizeInLimbs = 4, Index SuperblockSizeInLimbs = (1 << 16) / 64 / 4, Index SuperblockSelectStepSize = SuperblockSizeInLimbs * 64,
//         Index InSuperblockSelectStepSize = BlockSizeInLimbs * 64, typename SuperblockSelectT = U32, typename SuperblockRankT = Limb>
// class SelectBitvec
//     : public ClassicalBitvecImpl<SelectBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectStepSize, InSuperblockSelectStepSize, SuperblockSelectT, SuperblockRankT>,
//               BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectStepSize, InSuperblockSelectStepSize, SuperblockSelectT, SuperblockRankT> {

template<Index BlockSizeInLimbs = 4, Index SuperblockSizeInLimbs = (1 << 16) / 64 / 4, typename SuperblockSelectT = U32>
struct SelectBitvec
    : public ClassicalBitvecImpl<SelectBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectT>, BlockSizeInLimbs, SuperblockSizeInLimbs,
              SuperblockSizeInLimbs * 64, BlockSizeInLimbs * 64, SuperblockSelectT, Limb, OverwriteType<Byte>> {
    using Base = ClassicalBitvecImpl<SelectBitvec<BlockSizeInLimbs, SuperblockSizeInLimbs, SuperblockSelectT>, BlockSizeInLimbs,
            SuperblockSizeInLimbs, SuperblockSizeInLimbs * 64, BlockSizeInLimbs * 64, SuperblockSelectT, Limb, OverwriteType<Byte>>;
    using Base::Base;
    friend Base;
    static_assert(BlockSizeInLimbs <= 4, "BlockSize not supported");
    static_assert(BlockSizeInLimbs * 256 >= SuperblockSizeInLimbs, "SuperblockSize not supported");


    ADS_CPP20_CONSTEXPR void buildRankMetadata(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < this->numSuperblocks());
        if (superblockIdx == 0) {
            this->setSuperblockRank(0, 0);
        }
        Index inSuperblockSoFar = 0;
        Index beforeBlock = 0;
        Index numBlocks = this->numBlocksInSuperblock();
        if (superblockIdx == this->numSuperblocks() - 1) {
            numBlocks = this->numAccessibleBlocks() - numBlocks * superblockIdx;
        }
        Index superblockStartIdx = superblockIdx * this->numLimbsInSuperblock();
        for (Index b = 0; b < numBlocks; ++b) {
            this->setBlockRank_(superblockIdx, b, inSuperblockSoFar);
            for (Index i = 0; i < this->numLimbsInBlock(); ++i) {
                inSuperblockSoFar += popcount(this->getLimb(superblockStartIdx + this->numLimbsInBlock() * b + i));
            }
            if (inSuperblockSoFar / this->blockSize() > beforeBlock / this->blockSize()) {
            }
            beforeBlock = inSuperblockSoFar;
        }
        // there are numSuperblocks() + 1 super block counts
        this->setSuperblockRank(superblockIdx + 1, this->getSuperblockRank(superblockIdx) + inSuperblockSoFar);
    }


    constexpr SelectBitvec(UninitializedTag, Index numBits, CacheLine* mem = nullptr) noexcept
        : Base(UninitializedTag{}, numBits, mem) {}


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOneUnchecked(Index pos) const noexcept {
        return this->template rankFallback<true>(pos); // not fast, but that's not a problem
    }

    //    template<bool IsOne>
    //    [[nodiscard]] [[using gnu: hot]] ADS_CPP20_CONSTEXPR Index selectBlockIdx(Index& bitRank, Index superblockIdx) const noexcept {
    //        ADS_ASSUME(bitRank >= 0);
    //        ADS_ASSUME(bitRank < this->superblockSize());
    //        ADS_ASSUME(superblockIdx >= 0);
    //        ADS_ASSUME(superblockIdx < this->numSuperblocks());
    //        auto selectEntry = [](Index superblockIdx, Index inSuperblock) {
    //            ADS_ASSUME(inSuperblock >= 0);
    //            ADS_ASSUME(inSuperblock <= SelectBitvec::numBlocksInSuperblock());
    //            if constexpr (IsOne) {
    //                return superblockIdx * (SelectBitvec::numBlocksInSuperblock() + 2) + inSuperblock;
    //            } else {
    //                Index endOfBlock = (superblockIdx + 1) * (SelectBitvec::numBlocksInSuperblock() + 2) - 1;
    //                return endOfBlock - inSuperblock;
    //            }
    //        };
    //        Index inSuperblock = bitRank / this->blockSize();
    //        Index lowerSelectIdx = selectEntry(superblockIdx, inSuperblock);
    //        Index upperSelectIdx = selectEntry(superblockIdx, inSuperblock + 1);
    //        ADS_ASSUME(lowerSelectIdx >= 0);
    //        Index blockLowerBound = this->selectInSuperblock[lowerSelectIdx];
    //        Index blockUpperBound = this->selectInSuperblock[upperSelectIdx] + 1;
    //        ADS_ASSUME(blockLowerBound >= 0);
    //        ADS_ASSUME(blockUpperBound > blockLowerBound);
    //        ADS_ASSUME(blockUpperBound <= this->numBlocksInSuperblock());
    //        //        Index inSuperblockSoFar = IsOne ? lowerSelectIdx : SelectBitvec::numBlocksInSuperblock() + 1 - lowerSelectIdx;
    //        blockLowerBound += superblockIdx * this->numBlocksInSuperblock();
    //        blockUpperBound += superblockIdx * this->numBlocksInSuperblock();
    //        return this->template selectBlockIdxInRange<IsOne>(bitRank, superblockIdx, blockLowerBound, blockUpperBound);
    //    }
};

} // namespace ads

#endif // ADS_SELECT_BITVEC_HPP
