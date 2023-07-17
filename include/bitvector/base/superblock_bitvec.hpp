#ifndef ADS_SUPERBLOCK_BITVEC_HPP
#define ADS_SUPERBLOCK_BITVEC_HPP

#include "normal_bitvec.hpp"

namespace ads {

/// \brief CRTP base class of all bitvectors with fixed-sized rank blocks and superblocks. See BitvecBase for a CRTP
/// discussion. Derived classes should declare both this class and its base classes, available under the Base and
/// Base::Base alias, as friends. \tparam Derived The actual bitvector, which inherits from SuperblockBitvecBase<Derived>.
template<typename Derived, Index NumLimbsInSuperblock, Index NumLimbsInBlock, typename SuperblockRankT = Limb,
        Index NumLimbsInCacheLine = U64_PER_CACHELINE, typename OverwriteBlockRankT = FalseT>
class [[nodiscard]] SuperblockBitvecBase : public NormalBitvecBase<Derived, NumLimbsInBlock, NumLimbsInCacheLine> {
protected:
    friend Derived;

    using Base = NormalBitvecBase<Derived, NumLimbsInBlock, NumLimbsInCacheLine>;
    using Base::Base;
    using Base::derived;

    static_assert(NumLimbsInSuperblock > 0);
    static_assert(NumLimbsInSuperblock % NumLimbsInCacheLine == 0);
    static_assert(NumLimbsInCacheLine % NumLimbsInBlock == 0);

    // ** These aliases can be used unqualified (ie without Derived::) because they should refer to the actual types **
    using BlockRank
            = std::conditional_t<OverwriteBlockRankT{}, Unwrap<OverwriteBlockRankT>, IntType<bytesNeededForIndexing(NumLimbsInSuperblock * 64)>>;
    using SuperblockRank = SuperblockRankT;


public:
    // The following methods must be implemented by a base class
    //    [[nodiscard]] ADS_CPP20_CONSTEXPR BlockRank getBlockRank(Index i) const noexcept {}
    //
    //    [[nodiscard]] ADS_CPP20_CONSTEXPR SuperblockRank getSuperblockRank(Index i) const noexcept {}
    //
    //    ADS_CPP20_CONSTEXPR void setBlockRank(Index i, BlockRank newVal) noexcept {}
    //
    //    ADS_CPP20_CONSTEXPR void setBlockRank_(Index superblockIdx, Index blockInSuperblock, BlockRank newVal) noexcept {}
    //
    //    ADS_CPP20_CONSTEXPR void setSuperblockRank(Index i, SuperblockRank newVal) noexcept {}

public:
    // *** Default implementations for (static) getters ***

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numSuperblocks() const noexcept {
        return Derived::numSuperblocksForBits(derived().size());
    }

    [[nodiscard]] static constexpr Index superblockSize() noexcept { return NumLimbsInSuperblock * 64; }

    [[nodiscard]] static constexpr Index numLimbsInSuperblock() noexcept { return NumLimbsInSuperblock; }

    [[nodiscard]] static constexpr Index numCacheLinesInSuperblock() noexcept {
        ADS_ASSUME(Derived::numLimbsInSuperblock() % Derived::numLimbsInCacheLine() == 0);
        return Derived::numLimbsInSuperblock() / Derived::numLimbsInCacheLine();
    }

    [[nodiscard]] static constexpr Index numBlocksInSuperblock() noexcept {
        return Derived::numLimbsInSuperblock() / Derived::numLimbsInBlock();
    }

    [[nodiscard]] static constexpr Index bytesPerBlockRank() noexcept { return roundUpLog2(Derived::blockSize()); }

    [[nodiscard]] static constexpr Index bytesPerSuperblockRank() noexcept {
        return roundUpLog2(U64(Derived::superblockSize()));
    }

    [[nodiscard]] static constexpr Index numSuperblocksForBits(Index numBits) noexcept {
        return roundUpDiv(numBits, Derived::superblockSize());
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numOnes() const noexcept {
        return derived().getSuperblockRank(derived().numSuperblocks());
    }


    // *** Default implementations for Bitvector operations that are more than just (static) getters/setters ***

    // ** low-level operations that are nevertheless public because sometimes it's convenient or
    // faster to let the user call them **

    ADS_CPP20_CONSTEXPR void buildRankMetadata(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < derived().numSuperblocks());
        if (superblockIdx == 0) {
            derived().setSuperblockRank(0, 0);
        }
        Index inSuperblockSoFar = 0;
        Index numBlocks = derived().numBlocksInSuperblock();
        if (superblockIdx == derived().numSuperblocks() - 1) {
            numBlocks = derived().numAccessibleBlocks() - numBlocks * superblockIdx;
        }
        Index superblockStartIdx = superblockIdx * numLimbsInSuperblock();
        for (Index b = 0; b < numBlocks; ++b) {
            derived().setBlockRank_(superblockIdx, b, inSuperblockSoFar);
            for (Index i = 0; i < derived().numLimbsInBlock(); ++i) {
                inSuperblockSoFar += popcount(derived().getLimb(superblockStartIdx + derived().numLimbsInBlock() * b + i));
            }
        }
        // there are numSuperblocks() + 1 super block counts
        derived().setSuperblockRank(superblockIdx + 1, derived().getSuperblockRank(superblockIdx) + inSuperblockSoFar);
    }

    ADS_CPP20_CONSTEXPR void buildRankMetadata() noexcept {
        for (Index i = 0; i < derived().numSuperblocks(); ++i) {
            derived().buildRankMetadata(i);
        }
    }

    ADS_CPP20_CONSTEXPR void buildSelectMetadata() noexcept {
        // nothing to do in this default implementation
    }

    ADS_CPP20_CONSTEXPR void buildMetadata() noexcept {
        derived().completeWithZeros();
        derived().buildRankMetadata();
        derived().buildSelectMetadata();
    }

    /// Unlike rankOne(), this doesn't check that `pos` is valid, which gives no measurable performance benefits
    /// but makes implementations in subclasses slightly simpler as they don't need to perform bounds checking again.
    /// However, the combined ASSUME macros do improve performance by quite a bit (if the compiler couldn't assume that pos >= 0,
    //// performance would actually be significantly lower than with the throwing checks)
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOneUnchecked(Index pos) const noexcept {
        ADS_ASSUME(0 <= pos);
        ADS_ASSUME(pos < derived().numBits());
        ADS_ASSUME(pos >= 0);
        Index limbIdx = pos / 64;
        Index superblockIdx = limbIdx / derived().numLimbsInSuperblock();
        ADS_ASSUME(superblockIdx >= 0);
        Index blockIdx = limbIdx / derived().numLimbsInBlock();
        ADS_ASSUME(blockIdx >= 0);
        Index res = derived().getSuperblockRank(superblockIdx) + derived().getBlockRank(blockIdx);
        ADS_ASSUME(res >= 0);
        for (Index i = blockIdx * derived().numLimbsInBlock(); i < limbIdx; ++i) {
            ADS_ASSUME(limbIdx - i < derived().numLimbsInBlock());
            res += popcount(derived().getLimb(i));
        }
        return res + popcountBefore(derived().getLimb(limbIdx), pos % 64);
    }


protected:
    template<bool IsOne>
    [[nodiscard]] constexpr Index superblockRank(Index superblockIdx) const noexcept {
        ADS_ASSUME(superblockIdx >= 0);
        ADS_ASSUME(superblockIdx <= derived().numSuperblocks());
        Index rankOne = derived().getSuperblockRank(superblockIdx);
        ADS_ASSUME(rankOne >= 0);
        if constexpr (IsOne) {
            return rankOne;
        } else {
            // if superblockIdx is derived().numSuperblocks(), this can be greater than the number of zeros in the
            // bitvector, but that's not a problem
            Index numBitsBefore = superblockIdx * derived().superblockSize();
            ADS_ASSUME(numBitsBefore >= 0);
            ADS_ASSUME(rankOne <= numBitsBefore);
            return numBitsBefore - rankOne;
        }
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectSuperblockIdxInRange(Index& bitRank, Index l, Index u) const noexcept {
        ADS_ASSUME(l >= 0);
        ADS_ASSUME(l <= u);
        ADS_ASSUME(u > 0);
        ADS_ASSUME(u <= numSuperblocks());

        constexpr Index linearFallbackSize = 8;
        while (u - l > linearFallbackSize) [[unlikely]] {
            ADS_ASSUME(0 <= l);
            ADS_ASSUME(l <= u);
            Index mid = (l + u) / 2;
            Index midRank = superblockRank<IsOne>(mid);
            ADS_ASSUME(midRank >= superblockRank<IsOne>(l));
            if (midRank <= bitRank) {
                l = mid;
            } else {
                u = mid;
            }
        }
        ADS_ASSUME(u - l <= linearFallbackSize);
        ADS_ASSUME(l <= u);
        for (Index i = l; i < u; ++i) {
            if (superblockRank<IsOne>(i) > bitRank) {
                bitRank -= superblockRank<IsOne>(i - 1);
                return i - 1;
            }
        }
        ADS_ASSUME(superblockRank<IsOne>(u) >= bitRank);
        bitRank -= superblockRank<IsOne>(u - 1);
        return u - 1;
    }

    // ** internal implementations for the default select. They simply do binary searches, but the
    // constant factor is smaller than for inefficientSelect(). Still, some Bitvector implementations
    // provide faster select operations **
    template<bool IsOne>
    [[nodiscard]] constexpr Index selectSuperblockIdx(Index& bitRank) const noexcept {
        constexpr Index linearFallbackSize = 8;
        ADS_ASSUME(bitRank >= 0);
        // + 1 because we're searching for the first superblock where the rank is greater than bitRank
        Index l = bitRank / derived().superblockSize() + 1;
        Index u = derived().numSuperblocks();
        ADS_ASSUME(l > 0);
        ADS_ASSUME(l <= u);
        if (u - l > linearFallbackSize) [[unlikely]] {
            // set u close to the expected location for iid bit with 50% probability for '1', then increase exponentially
            // until it is an upper bound. Unlike binary search, this starts with a less pessimistic search window. This improves performance for random values but hurts for especially hard cases.
            u = l;
            do {
                l = u;
                u *= 2;
                if (u >= derived().numSuperblocks()) {
                    u = derived().numSuperblocks();
                    break;
                }
            } while (superblockRank<IsOne>(u) <= bitRank);
        }
        return derived().template selectSuperblockIdxInRange<IsOne>(bitRank, l, u);
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index blockRank(Index blockIdx) const noexcept {
        ADS_ASSUME(blockIdx >= 0);
        Index rankOne = derived().getBlockRank(blockIdx);
        ADS_ASSUME(rankOne >= 0);
        ADS_ASSUME(rankOne <= derived().numOnes());
        if constexpr (IsOne) {
            return rankOne;
        } else {
            Index numBitsBefore = (blockIdx % derived().numBlocksInSuperblock()) * derived().blockSize();
            ADS_ASSUME(rankOne <= numBitsBefore);
            return numBitsBefore - rankOne;
        }
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectBlockIdxInRange(Index& bitRank, Index superblockIdx, Index l, Index u) const noexcept {
        constexpr Index linearFallbackSize = 2 * sizeof(Limb) / sizeof(BlockRank);
        ADS_ASSUME(u >= l);
        ADS_ASSUME(u - l <= derived().numBlocksInSuperblock());
        ADS_ASSUME(u / this->numBlocksInSuperblock() - l / this->numBlocksInSuperblock() <= 1);
        while (u - l > linearFallbackSize) {
            Index mid = (l + u) / 2;
            Index midRank = blockRank<IsOne>(mid);
            if (midRank > bitRank) {
                u = mid;
            } else {
                l = mid;
            }
        }
        ADS_ASSUME(l >= 0);
        ADS_ASSUME(u >= l);
        ADS_ASSUME(u - l <= linearFallbackSize);
        ADS_ASSUME(l >= superblockIdx * derived().numBlocksInSuperblock());
        for (Index i = l; i < u; ++i) { // TODO: Start from l + 1
            ADS_ASSUME((i % derived().numBlocksInSuperblock()) == 0 || blockRank<IsOne>(i) >= blockRank<IsOne>(i - 1));
            ADS_ASSUME((i % derived().numBlocksInSuperblock()) == 0 || blockRank<IsOne>(i - 1) <= bitRank);
            if (blockRank<IsOne>(i) > bitRank) {
                bitRank -= blockRank<IsOne>(i - 1);
                return i - 1;
            }
        }
        bitRank -= blockRank<IsOne>(u - 1);
        ADS_ASSUME(bitRank >= 0);
        ADS_ASSUME(bitRank < this->blockSize());
        return u - 1;
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectBlockIdx(Index& bitRank, Index superblockIdx) const noexcept {
        ADS_ASSUME(superblockIdx >= 0);
        ADS_ASSUME(superblockIdx < derived().numSuperblocks());
        ADS_ASSUME(bitRank >= 0);
        ADS_ASSUME(bitRank < derived().superblockSize());
        // we're searching for the first block with count strictly greater than bitRank
        Index l = superblockIdx * derived().numBlocksInSuperblock() + 1;
        Index u = std::min(l + derived().numBlocksInSuperblock() - 1, derived().numBlocks());
        return derived().template selectBlockIdxInRange<IsOne>(bitRank, superblockIdx, l, u);
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectBitIdx(Limb limb, Index bitIndex) const noexcept {
        if constexpr (IsOne) {
            return u64Select(limb, bitIndex);
        } else {
            return u64Select(~limb, bitIndex);
        }
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectInBlock(Index bitRank, Index blockIdx) const noexcept {
        if constexpr (NumLimbsInBlock == 1) {
            Index limb = this->getLimb(blockIdx);
            return selectBitIdx<IsOne>(limb, bitRank);
        } else if constexpr (NumLimbsInBlock == 4) {
            if constexpr (IsOne) {
                return u256Select(this->getBlock(blockIdx).data(), bitRank);
            } else {
                return u256SelectZero(this->getBlock(blockIdx).data(), bitRank);
            }
        } else {
            auto rankFunc = [](Limb limb) noexcept {
                if constexpr (IsOne) {
                    return popcount(limb);
                } else {
                    return 64 - popcount(limb);
                }
            };
            Span block = this->getBlock(blockIdx);
            for (Index i = 0; i + 1 < NumLimbsInBlock; ++i) {
                Limb limb = block[i];
                Index rank = rankFunc(limb);
                ADS_ASSUME(rank >= 0);
                if (rank > bitRank) {
                    return i * 64 + selectBitIdx<IsOne>(limb, bitRank);
                }
                bitRank -= rank;
                ADS_ASSUME(bitRank >= 0);
            }
            Limb limb = block[NumLimbsInBlock - 1];
            ADS_ASSUME(bitRank >= 0);
            ADS_ASSUME(rankFunc(limb) > bitRank);
            return (NumLimbsInBlock - 1) * 64 + selectBitIdx<IsOne>(limb, bitRank);
        }
    }

public:
    template<bool IsOne>
    [[nodiscard]] constexpr Index select(Index bitRank) const {
        if (bitRank < 0 || bitRank >= derived().size()) [[unlikely]] {
            ADS_THROW("invalid rank for select query: " + std::to_string(bitRank));
        }
        Index superblockIdx = derived().template selectSuperblockIdx<IsOne>(bitRank);
        ADS_ASSUME(superblockIdx >= 0);
        ADS_ASSUME(superblockIdx < derived().numSuperblocks());
        ADS_ASSUME(bitRank >= 0);
        ADS_ASSUME(bitRank < this->superblockSize());
        Index blockIdx = derived().template selectBlockIdx<IsOne>(bitRank, superblockIdx);
        ADS_ASSUME(blockIdx >= 0);
        ADS_ASSUME(blockIdx < derived().numBlocks());
        ADS_ASSUME(bitRank >= 0);
        ADS_ASSUME(bitRank < this->blockSize());
        Index inBlock = derived().template selectInBlock<IsOne>(bitRank, blockIdx);
        ADS_ASSUME(inBlock >= 0);
        ADS_ASSUME(inBlock < NumLimbsInBlock * 64);
        return blockIdx * derived().blockSize() + inBlock;
    }

    //    [[nodiscard]] constexpr std::pair<Index, Index> selectOneAndPrevOne(Index secondRank) const {
    //        if (secondRank < 0 || secondRank >= derived().size()) [[unlikely]] {
    //            ADS_THROW("invalid rank for select query: " + std::to_string(secondRank));
    //        }
    //        Index superblockIdx = derived().template selectSuperblockIdx<true>(secondRank);
    //        ADS_ASSUME(superblockIdx >= 0);
    //        ADS_ASSUME(superblockIdx < derived().numSuperblocks());
    //        Index blockIdx = derived().template selectBlockIdx<true>(secondRank, superblockIdx);
    //        ADS_ASSUME(blockIdx >= 0);
    //        ADS_ASSUME(blockIdx < derived().numBlocks());
    //        Index inBlock = derived().template selectInBlock<true>(secondRank, blockIdx);
    //        ADS_ASSUME(inBlock >= 0);
    //        ADS_ASSUME(inBlock < NumLimbsInBlock * 64);
    //        return blockIdx * derived().blockSize() + inBlock;
    //    }
};


} // namespace ads

#endif // ADS_SUPERBLOCK_BITVEC_HPP
