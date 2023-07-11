#ifndef ADS_SUPERBLOCK_BITVEC_HPP
#define ADS_SUPERBLOCK_BITVEC_HPP

#include "normal_bitvec.hpp"

namespace ads {

/// \brief CRTP base class of all bitvectors with fixed-sized rank blocks and superblocks. See BitvecBase for a CRTP
/// discussion. Derived classes should declare both this class and its base classes, available under the Base and
/// Base::Base alias, as friends. \tparam Derived The actual bitvector, which inherits from SuperblockBitvecBase<Derived>.
template<typename Derived, Index NumLimbsInSuperblock, Index NumLimbsInBlock, typename SuperblockRankT = Limb, Index NumLimbsInCacheLine = U64_PER_CACHELINE>
class [[nodiscard]] SuperblockBitvecBase : public NormalBitvecBase<Derived, NumLimbsInBlock, NumLimbsInCacheLine> {
    friend Derived;

    using Base = NormalBitvecBase<Derived, NumLimbsInBlock, NumLimbsInCacheLine>;
    using Base::Base;
    using Base::derived;

    static_assert(NumLimbsInSuperblock > 0);
    static_assert(NumLimbsInSuperblock % NumLimbsInCacheLine == 0);
    static_assert(NumLimbsInCacheLine % NumLimbsInBlock == 0);

    // ** These aliases can be used unqualified (ie without Derived::) because they should refer to the actual types **
    using BlockRank = IntType<bytesNeededForIndexing(NumLimbsInSuperblock * 64)>;
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
            numBlocks = derived().numBlocks() - numBlocks * superblockIdx;
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
        if (pos <= 0) [[unlikely]] {
            return 0;
        }
        --pos;
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
        return res + popcountUntil(derived().getLimb(limbIdx), pos % 64);
    }


private:
    // ** internal implementations for the default select. They simply do binary searches, but the
    // constant factor is smaller than for inefficientSelect(). Still, some Bitvector implementations
    // provide faster select operations **
    template<bool IsOne>
    [[nodiscard]] constexpr Index selectSuperBlockIdx(Index& bitRank) const {
        constexpr Index linearFallbackSize = 8;
        auto rankFunc = [this](Index i) noexcept {
            ADS_ASSUME(i >= 0);
            ADS_ASSUME(i <= derived().numSuperblocks());
            Index rankOne = derived().getSuperblockRank(i);
            ADS_ASSUME(rankOne >= 0);
            if constexpr (IsOne) {
                return rankOne;
            } else {
                // if i is derived().numSuperblocks(), this can be greater than the number of zeros in the
                // bitvector, but that's not a problem
                Index numBitsBefore = i * derived().superblockSize();
                ADS_ASSUME(rankOne <= numBitsBefore);
                return numBitsBefore - rankOne;
            }
        };
        ADS_ASSUME(bitRank >= 0);
        // + 1 because we're searching for the first superblock where the rank is greater than bitRank
        Index l = bitRank / derived().superblockSize() + 1;
        Index u = derived().numSuperblocks();
        ADS_ASSUME(l > 0);
        ADS_ASSUME(l <= u);
        if (u - l > linearFallbackSize) {
            // set u close to the expected location for iid bit with 50% probability for '1', then increase exponentially
            // until it is an upper bound. Unlike binary search, this starts with a less pessimistic search window and
            // should hopefully be easier on the branch predictor. This improves performance for random values but hurts for especially hard cases.
            u = l;
            do {
                l = u;
                u *= 2;
                if (u >= derived().numSuperblocks()) {
                    u = derived().numSuperblocks();
                    break;
                }
            } while (rankFunc(u) <= bitRank);
            while (u - l > linearFallbackSize) {
                ADS_ASSUME(0 < l);
                ADS_ASSUME(l <= u);
                Index mid = (l + u) / 2;
                Index midRank = rankFunc(mid);
                ADS_ASSUME(midRank >= rankFunc(l));
                if (midRank <= bitRank) {
                    l = mid;
                } else {
                    u = mid;
                }
            }
        }
        ADS_ASSUME(u - l <= linearFallbackSize);
        ADS_ASSUME(l <= u);
        for (Index i = l; i < u; ++i) {
            if (rankFunc(i) > bitRank) {
                bitRank -= rankFunc(i - 1);
                return i - 1;
            }
        }
        ADS_ASSUME(rankFunc(u) >= bitRank);
        bitRank -= rankFunc(u - 1);
        return u - 1;
    }

    // TODO: For most bitvectors, superblock indices can be represented with 32 bit values.
    template<bool IsOne>
    [[nodiscard]] constexpr Index selectBlockIdx(Index& bitRank, Index superBlockIdx) const noexcept {
        auto rankFunc = [this](Index i) noexcept {
            Index rankOne = derived().getBlockRank(i);
            if constexpr (IsOne) {
                return rankOne;
            } else {
                Index numBitsBefore = (i % derived().numBlocksInSuperblock()) * derived().blockSize();
                ADS_ASSUME(rankOne <= numBitsBefore);
                return numBitsBefore - rankOne;
            }
        };
        constexpr Index linearFallbackSize = 2 * sizeof(Limb) / sizeof(BlockRank);
        ADS_ASSUME(superBlockIdx >= 0);
        ADS_ASSUME(superBlockIdx < derived().numSuperblocks());
        ADS_ASSUME(bitRank >= 0);
        // we're searching for the first block with count strictly greater than bitRank
        Index l = superBlockIdx * derived().numBlocksInSuperblock() + 1;
        Index u = std::min(l + derived().numBlocksInSuperblock() - 1, derived().numBlocks());
        ADS_ASSUME(u >= l);
        ADS_ASSUME(u - l < derived().numBlocksInSuperblock());
        while (u - l > linearFallbackSize) {
            Index mid = (l + u) / 2;
            Index midRank = rankFunc(mid);
            if (midRank > bitRank) {
                u = mid;
            } else {
                l = mid;
            }
        }
        ADS_ASSUME(u >= l);
        ADS_ASSUME(u - l <= linearFallbackSize);
        ADS_ASSUME(l > superBlockIdx * derived().numBlocksInSuperblock());
        for (Index i = l; i < u; ++i) {
            ADS_ASSUME(i == 0 || derived().getBlockRank(i) >= derived().getBlockRank(i - 1));
            if (rankFunc(i) > bitRank) {
                bitRank -= rankFunc(i - 1);
                return i - 1;
            }
        }
        bitRank -= rankFunc(u - 1);
        return u - 1;
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectLimbIdx(Index& bitRank, Index blockIdx) const noexcept {
        Index first = blockIdx * derived().numLimbsInBlock();
        auto rankFunc = [this](Index i) noexcept {
            Limb l = derived().getLimb(i);
            if constexpr (IsOne) {
                return popcount(l);
            } else {
                return 64 - popcount(l);
            }
        };
        for (Index i = first; i < derived().numLimbs(); ++i) {
            Index rank = rankFunc(i);
            ADS_ASSUME(i - first < derived().numLimbsInBlock());
            ADS_ASSUME(rank >= 0);
            if (rank > bitRank) {
                return i;
            }
            bitRank -= rank;
            ADS_ASSUME(bitRank >= 0);
        }
        return derived().numLimbs() - 1;
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectBitIdx(Limb limb, Index bitIndex) const noexcept {
        if constexpr (IsOne) {
            return u64Select(limb, bitIndex);
        } else {
            return u64Select(~limb, bitIndex);
        }
    }

public:
    template<bool IsOne>
    [[nodiscard]] constexpr Index select(Index bitRank) const {
        if (bitRank < 0 || bitRank >= derived().size()) [[unlikely]] {
            throw std::invalid_argument("invalid rank for select query: " + std::to_string(bitRank));
        }
        Index superBlockIdx = derived().template selectSuperBlockIdx<IsOne>(bitRank);
        Index blockIdx = derived().template selectBlockIdx<IsOne>(bitRank, superBlockIdx);
        Index limbIdx = derived().template selectLimbIdx<IsOne>(bitRank, blockIdx);
        Index bitIdx = derived().template selectBitIdx<IsOne>(derived().getLimb(limbIdx), bitRank);
        return limbIdx * 64 + bitIdx;
    }
};


} // namespace ads

#endif // ADS_SUPERBLOCK_BITVEC_HPP
