#ifndef ADS_RECURSIVE_BITVEC_HPP
#define ADS_RECURSIVE_BITVEC_HPP

#include "bitvec_base.hpp"
#include "trivial_bitvec.hpp"

namespace ads {

// ZERO_ONLY isn't implemented as it can be emulated by storing every bit flipped
enum class SupportedSelects { BOTH, ONE_ONLY };

/// \brief Recursive Bitvector. Most of the later template arguments shouldn't really need to be changed.
/// Needs about 4% space overhead (compared to around 3% space overhead for the EfficientRankBitvec) but answers select
/// queries faster starting at about 10^5 bits.
/// \tparam NestedBitvec The nested bitvector, of which two with approx. size() / BlockSize bits are stored.
/// \tparam SelectOps Which select operations to support. Not supporting both saves one of two bitvectors, which
/// stores approx. size()/BlockSize logical bits (its allocated size may be larger, depending on the type NestedBitvec).
/// \tparam BlockSize The number of bits in a block. This is used to determines the BlockCountT type;
/// smaller values use more space but can reduce query times. Must currently be a multiple of 64.
/// \tparam BlockCountT The type used to represent block counts, must be able to hold values up to BlockSize - 1.
/// Shouldn't really need to be overwritten.
template<ADS_BITVEC_CONCEPT NestedBitvec, SupportedSelects SelectOps = SupportedSelects::BOTH, Index BlockSize = 256,
        typename BlockCountT = IntType<bytesNeededForIndexing(BlockSize)>>
class RecursiveBitvec : public NormalBitvecBase<RecursiveBitvec<NestedBitvec, SelectOps, BlockSize, BlockCountT>> {
    static_assert(sizeof(BlockCountT) >= bytesNeededForIndexing(BlockSize));
    static_assert(BlockSize >= 64);

    using Base = NormalBitvecBase<RecursiveBitvec<NestedBitvec, SelectOps, BlockSize, BlockCountT>>;
    friend Base;
    friend Base::Base;
    using Nested = NestedBitvec;

    /// \brief The nestedZeros Bitvector stores a bit for each block b, which is one iff the block contains a zero at
    /// index i in [b * BlockSize, (b+1) * BlockSize) with rankZero(i) % BlockSize == 0 and rankZero(i) > 0.
    /// The nestedOnes Bitvector does the same for ones.
    Nested nestedZeros = Nested(); // store zeros bitvector before ones bitvector to keep the ones bitvector closer
    Nested nestedOnes = Nested();  // to the blockRanks array for better data locality on rank / selectOne
    /// \brief The blockRanks array stores the number of 1 at the end of the respective block, modulo BlockSize.
    View<BlockCountT> blockRanks = View<BlockCountT>();
    /// \brief The actual bit sequence.
    View<Limb> vec = View<Limb>();
    /// \brief The logical number of bits in the bit sequence. The bitvector may internally allocate, set and access
    /// a greater number of bits that this value.
    Index numBits = 0;

    constexpr static Index limbsInBlock = BlockSize / 64;

    constexpr static Index nestedSizeForBits(Index numBits) noexcept {
        // append a "1" to the nested bitvector to avoid a special case in select
        // also, the extra entry in the blockRanks array is additionally used for numOnes()
        return roundUpDiv(numBits, BlockSize) + 1;
    }

    ADS_CPP20_CONSTEXPR RecursiveBitvec(UninitializedTag, Index numBits, Limb* mem)
        : Base(numBits, mem), numBits(numBits) {
        Index nestedSize = nestedSizeForBits(numBits); // for the blockRanks array, the last entry can be used for numOnes()
        Index offset = 0;
        if constexpr (SelectOps != SupportedSelects::ONE_ONLY) {
            nestedZeros = Nested::uninitializedForSize(nestedSize, this->allocation.memory());
            offset += nestedZeros.allocatedSizeInLimbs();
            ADS_ASSUME(nestedSize == nestedZeros.size());
        }
        nestedOnes = Nested::uninitializedForSize(nestedSize, this->allocation.memory() + offset);
        offset += nestedOnes.allocatedSizeInLimbs();
        ADS_ASSUME(nestedSize == nestedOnes.size());
        blockRanks = View<BlockCountT>(this->allocation.memory() + offset, nestedSize);
        Index numLimbsInBitSequence = roundUpDiv(numBits, BlockSize) * limbsInBlock; // round up to full blocks
        vec = View<Limb>(this->allocation.memory() + this->allocation.size() - numLimbsInBitSequence, numLimbsInBitSequence);
        ADS_ASSUME(offset + blockRanks.sizeInLimbs() + vec.sizeInLimbs() == this->allocation.size());
    }


public:
    constexpr RecursiveBitvec() noexcept = default;

    explicit constexpr RecursiveBitvec(Index numBits, Limb fill, Limb* mem = nullptr)
        : RecursiveBitvec(UninitializedTag{}, numBits, mem) {
        this->fill(fill);
    }

    explicit ADS_CPP20_CONSTEXPR RecursiveBitvec(Span<const Limb> limbs, Limb* mem = nullptr) noexcept
        : RecursiveBitvec(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR RecursiveBitvec(Span<const Limb> limbs, Index numBits, Limb* mem = nullptr) noexcept
        : RecursiveBitvec(UninitializedTag{}, numBits, mem) {
        this->copyFrom(limbs);
    }

    explicit ADS_CPP20_CONSTEXPR RecursiveBitvec(std::string_view str, Index base = 2, Limb* mem = nullptr) noexcept
        : RecursiveBitvec(UninitializedTag{}, str.size() * intLog2(base), mem) {
        this->initFromStr(str, base);
    }

    [[nodiscard]] static constexpr Index allocatedSizeInLimbsForBits(Index numBits) noexcept {
        Index numBlocks = nestedSizeForBits(numBits);
        const Index factor = (SelectOps == SupportedSelects::BOTH ? 2 : 1);
        Index nestedSizesInLimbs = Nested::allocatedSizeInLimbsForBits(numBlocks) * factor;
        Index arrSizeInLimbs = roundUpDiv(numBlocks * sizeof(BlockCountT), sizeof(Limb));
        Index bitSequenceSizeInLimbs = roundUpDiv(numBits, BlockSize) * limbsInBlock; // no incomplete blocks
        return nestedSizesInLimbs + arrSizeInLimbs + bitSequenceSizeInLimbs;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return numBits; }

    [[nodiscard]] constexpr Index numLimbs() const noexcept {
        ADS_ASSUME(vec.numT >= roundUpDiv(numBits, 64));
        return roundUpDiv(numBits, 64); // vec.numT rounds up towards full blocks
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBlocks() const noexcept { return nestedOnes.size() - 1; }

private:
    [[nodiscard]] ADS_CPP20_CONSTEXPR const Limb* getLimbArray() const noexcept { return vec.ptr; }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb* getLimbArray() noexcept { return vec.ptr; }

public:
    ADS_CPP20_CONSTEXPR void buildMetadata() noexcept {
        [[maybe_unused]] Index rankOne = 0;
        [[maybe_unused]] Index rankZero = 0;
        Index numBlocks = this->numBlocks();
        assert(numBlocks == roundUpDiv(numLimbs(), limbsInBlock));
        if (numLimbs() > 0) [[likely]] {
            Index numIgnoredBits = (64 - (size() % 64)) % 64;
            this->setLimb(numLimbs() - 1, (this->getLimb(numLimbs() - 1) << numIgnoredBits) >> numIgnoredBits);
        }
        for (Index i = numLimbs(); i < numBlocks * limbsInBlock; ++i) {
            vec.setBits(i, 0); // can't use setLimb() because that assumed that i < numLimbs()
        }

        for (Index blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
            Index prevRankOne = rankOne;
            Limb* blockPtr = &this->getLimbRef(blockIdx * limbsInBlock);
            rankOne += rankInBlock(blockPtr);
            nestedOnes.setBit(blockIdx, rankOne / BlockSize > prevRankOne / BlockSize);
            blockRanks.ptr[blockIdx] = rankOne % BlockSize;
            if constexpr (SelectOps != SupportedSelects::ONE_ONLY) {
                Index prevRankZero = rankZero;
                rankZero = (blockIdx + 1) * BlockSize - rankOne;
                nestedZeros.setBit(blockIdx, rankZero / BlockSize > prevRankZero / BlockSize);
            }
        }
        ADS_ASSUME(rankOne <= size());
        blockRanks.ptr[numBlocks] = rankOne % BlockSize; // makes numOnes() simpler and avoids special cases for select(largeVal)
        nestedOnes.setBit(numBlocks);
        nestedOnes.buildMetadata();
        if constexpr (SelectOps != SupportedSelects::ONE_ONLY) {
            nestedZeros.setBit(numBlocks); // more special case avoidance
            nestedZeros.buildMetadata();
        }
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numOnes() const noexcept {
        // The nested bitvector stores an additional 1;
        return (nestedOnes.numOnes() - 1) * BlockSize + blockRanks[numBlocks()];
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOneUnchecked(Index idx) const noexcept {
        ADS_ASSUME(idx >= 0);
        Index blockIdx = idx / BlockSize;
        ADS_ASSUME(blockIdx >= 0);
        ADS_ASSUME(blockIdx < numBlocks());
        Index nestedRank = nestedOnes.rankOne(blockIdx + 1); // + 1 because rankOne doesn't consider the current bit
        Index rankAfterBlock = blockRanks[blockIdx] + BlockSize * nestedRank;
        Index limbIdx = idx / 64;
        ADS_ASSUME(limbIdx >= blockIdx * limbsInBlock);
        ADS_ASSUME(limbIdx < limbsInBlock * (blockIdx + 1));
        Index res = rankAfterBlock;
        Index blockEndIdx = std::min(numLimbs(), (blockIdx + 1) * limbsInBlock);
        // count the number of 1s between idx and the end of the block. This doesn't require looking at
        // the previous block rank, which helps a lot when blockIdx is 0.
        for (Index i = blockEndIdx - 1; i > limbIdx; --i) { // TODO: Unroll? SIMD? Implement same for rankZero?
            res -= popcount(this->getLimb(i));
        }
        const Limb mask = (Limb(-1) << (idx % 64));
        return res - popcount(this->getLimb(limbIdx) & mask);
    }

    // Don't iterate over the blockRanks twice and don't call nested selectOne twice unless necessary.
    // Cf. rankOneUnchecked() below. This function exists because it fits the pattern in which Elias-Fano calls selectOne
    // and the pattern in which this bitvector calls the nested bitvector's selectOne function, so for a
    // RecursiveBitvec<RecursiveBitvec<TrivialBitvec>>>, this turns 14 selectOne()s per predecessor into usually
    // no more than 3 selectOneAndPrev calls overall. // TODO: Make EliasFano use selectOneAndPrevOne instead of selectZero.
    [[nodiscard]] ADS_CPP20_CONSTEXPR std::pair<Index, Index> selectOneAndPrevOne(Index rankOfSecond) const noexcept {
        ADS_ASSUME(rankOfSecond >= 0);
        ADS_ASSUME(rankOfSecond < this->numOnes());
        if (rankOfSecond == 0) [[unlikely]] {
            return {0, selectOneUnchecked(rankOfSecond)};
        }
        Index nestedSearchRank = rankOfSecond / BlockSize;
        if ((rankOfSecond - 1) / BlockSize != nestedSearchRank) [[unlikely]] {
            return {selectOneUnchecked(rankOfSecond - 1), selectOneUnchecked(rankOfSecond)};
        }
        Index upperRankInBlock = rankOfSecond % BlockSize;
        Index lowerRankInBlock = (rankOfSecond - 1) % BlockSize;
        auto [nestedIdx, nestedEndIdx] = nestedOnes.selectOneAndPrevOne(nestedSearchRank);
        Index blockIdx = selectBlock<true>(nestedIdx, nestedEndIdx, lowerRankInBlock);
        ADS_ASSUME(blockRanks[blockIdx] >= upperRankInBlock || blockIdx == nestedEndIdx);
        if (blockRanks[blockIdx] == upperRankInBlock) [[unlikely]] {
            // the second one is in a later block *this was the last one in the block), so do binary search again to find it
            Index upperBlock = selectBlock<true>(blockIdx + 1, nestedEndIdx, upperRankInBlock);
            Index rankInUpperBlockFromRight = (blockRanks[upperBlock] + BlockSize - upperRankInBlock - 1) % BlockSize;
            ADS_ASSUME(rankInUpperBlockFromRight >= 0);
            Index upperRes = selectInBlockFromRight<true>(upperBlock, rankInUpperBlockFromRight);

            [[maybe_unused]] Index rankInLowerBlockFromRight = (blockRanks[blockIdx] + BlockSize - lowerRankInBlock - 1) % BlockSize;
            ADS_ASSUME(rankInLowerBlockFromRight == 0);
            Index lowerRes = selectInBlockFromRight<true>(blockIdx, 0);
            return {lowerRes, upperRes};
        }
        ADS_ASSUME(blockIdx >= nestedIdx);

        Index rankInBlockFromRight = (blockRanks[blockIdx] + BlockSize - upperRankInBlock - 1) % BlockSize;
        ADS_ASSUME(rankInBlockFromRight >= 0);
        ADS_ASSUME(rankInBlockFromRight < rankInBlock(&this->getLimbRef(blockIdx * limbsInBlock)));
        return select2InBlockFromRight(blockIdx, rankInBlockFromRight);
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectOneUnchecked(Index rank) const noexcept {
        ADS_ASSUME(rank >= 0);
        ADS_ASSUME(rank < numOnes());
        Index nestedSearchRank = rank / BlockSize;
        rank %= BlockSize;
        ADS_ASSUME(nestedSearchRank < nestedOnes.numOnes());
        auto [nestedIdx, nestedEndIdx] = nestedOnes.selectOneAndPrevOne(nestedSearchRank);
        ADS_ASSUME(nestedIdx <= nestedEndIdx);
        Index blockIdx = selectBlock<true>(nestedIdx, nestedEndIdx, rank);
        ADS_ASSUME(blockIdx >= nestedIdx);
        // blockRanks[blockIdx] contains the number of 1s relative to nestedSearchRank * BlockSize or (nestedSearchRank +
        // 1) * BlockSize, so we need to find the (blockRanks[blockIdx] - rank mod BlockSize)th bit counting from the right, 0 indexed
        Index rankInBlockFromRight = (blockRanks[blockIdx] + BlockSize - rank - 1) % BlockSize;
        ADS_ASSUME(rankInBlockFromRight >= 0);
        ADS_ASSUME(rankInBlockFromRight < rankInBlock(&this->getLimbRef(blockIdx * limbsInBlock)));
        return selectInBlockFromRight<true>(blockIdx, rankInBlockFromRight);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectZeroUnchecked(Index rankZero) const noexcept {
        if constexpr (SelectOps == SupportedSelects::ONE_ONLY) {
            return this->template selectFallback<false>(rankZero);
        } else {
            return selectZeroUncheckedImpl(rankZero);
        }
    }

private:
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectZeroUncheckedImpl(Index rankZero) const noexcept {
        ADS_ASSUME(rankZero >= 0);
        Index nestedSearchRank = rankZero / BlockSize;
        rankZero %= BlockSize;
        ADS_ASSUME(nestedSearchRank < nestedZeros.numOnes());
        auto [nestedIdx, nestedEndIdx] = nestedZeros.selectOneAndPrevOne(nestedSearchRank);
        Index blockIdx = selectBlock<false>(nestedIdx, nestedEndIdx, rankZero);
        ADS_ASSUME(blockIdx >= nestedIdx);
        ADS_ASSUME(blockIdx < numBlocks());
        // blockRanks[blockIdx] contains the number of 1s relative to nestedSearchRan * BlockSize or
        // (nestedSearchRank + 1) * BlockSize, so we need to find the (blockRans[blockIdx] - rankZero mod
        // BlockSize)th bit counting from the right, 0 indexed
        Index rankInBlockFromRight = (2 * BlockSize - blockRanks[blockIdx] - rankZero - 1) % BlockSize;
        if (blockIdx == numBlocks() - 1) [[unlikely]] {
            // TODO: Store complete blocks, make sure all limbs after numLimbs() are 0
            Index numIgnoredZeros = ((limbsInBlock - numLimbs() % limbsInBlock) % limbsInBlock) * 64;
            rankInBlockFromRight -= numIgnoredZeros;
        }
        ADS_ASSUME(rankInBlockFromRight >= 0);
        ADS_ASSUME(rankInBlockFromRight < BlockSize);
        return selectInBlockFromRight<false>(blockIdx, rankInBlockFromRight);
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectInBlockFromRight(Index blockIdx, Index rank) const noexcept {
        ADS_ASSUME(blockIdx >= 0);
        ADS_ASSUME(rank >= 0);
        ADS_ASSUME(rank < BlockSize);
        Index blockStart = blockIdx * limbsInBlock;
        Index blockEnd = std::min((blockIdx + 1) * limbsInBlock, numLimbs());
        ADS_ASSUME(blockStart < blockEnd);
        Index rankInBlock = rank;
        auto limbValue = [this](Index i) {
            if constexpr (IsOne) {
                return this->getLimb(i);
            } else {
                return ~this->getLimb(i);
            }
        };
        for (Index i = blockEnd - 1; i >= blockStart; --i) {
            Limb limb = limbValue(i);
            Index popcnt = popcount(limb);
            if (rankInBlock < popcnt) {
                ADS_ASSUME(rankInBlock >= 0);
                Index bitIdxFromRight = popcnt - rankInBlock - 1;
                return i * 64 + u64Select(limb, bitIdxFromRight);
            }
            rankInBlock -= popcnt;
        }
        ADS_ASSUME(false);
        return -1;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR std::pair<Index, Index> select2InBlockFromRight(Index blockIdx, Index rankOfRight) const noexcept {
        ADS_ASSUME(blockIdx >= 0);
        ADS_ASSUME(rankOfRight >= 0);
        ADS_ASSUME(rankOfRight + 1 < BlockSize); // the left value's rank must also be in the same block
        Index blockStart = blockIdx * limbsInBlock;
        Index blockEnd = std::min((blockIdx + 1) * limbsInBlock, numLimbs());
        ADS_ASSUME(blockStart < blockEnd);
        Index rankInBlock = rankOfRight;
        for (Index i = blockEnd - 1; i >= blockStart; --i) {
            Limb limb = this->getLimb(i);
            Index popcnt = popcount(limb);
            if (rankInBlock < popcnt) {
                ADS_ASSUME(rankInBlock >= 0);
                Index rankInLimb = popcnt - rankInBlock - 1;
                // TODO: Use separate function for computing the rank of two values at the same time?
                // Should be relatively easy for the table lookup implementation.
                Index rightRes = i * 64 + u64Select(limb, rankInLimb);
                if (rankInBlock + 1 < popcnt) [[likely]] {
                    Index leftRes = i * 64 + u64Select(limb, rankInLimb - 1);
                    return {leftRes, rightRes};
                }
                while (this->getLimb(--i) == 0) {
                    ADS_ASSUME(i >= blockStart);
                }
                limb = this->getLimb(i);
                ADS_ASSUME(i >= blockStart);
                return {i * 64 + u64Select(limb, popcount(limb) - 1), rightRes};
            }
            rankInBlock -= popcnt;
        }
        ADS_ASSUME(false);
        return {-1, -1};
    }


    static constexpr Index linearFallbackSize = 16;

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index getRank(Index i) const noexcept {
        Index rank = blockRanks[i];
        ADS_ASSUME(rank >= 0);
        ADS_ASSUME(rank < BlockSize);
        if constexpr (IsOne) {
            return rank;
        } else {
            // The blockRank array stores rankOne() after the end of the block, modulo BlockSize.
            // rankZero((i + 1) * BlockSize) % BlockSize is ((i + 1) * BlockSize - blockRank[i]) % BlockSize, which
            // is just BlockSize - blockRank[i], except when the result would be BlockSize; it should be 0 instead.
            return (BlockSize - rank) % BlockSize;
        }
    }

    /// \brief The subrange [lower, upper) of block ranks must be in ascending order. Finds the first block in the range
    /// [lower, upper) where the rank is strictly larger than rank or returns upper if no such block exists.
    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectBlock(Index lower, Index upper, Index rank) const noexcept {
        ADS_ASSUME(rank >= 0);
        ADS_ASSUME(rank < BlockSize);
        ADS_ASSUME(lower >= 0);
        ADS_ASSUME(upper >= lower);
        ADS_ASSUME(upper <= numBlocks());
        //        // TODO: Remove this...
        //        auto compareRank = [](Index rankZero, BlockCountT blockRankOne) {
        //            if constexpr (IsOne) {
        //                return rankZero < blockRankOne;
        //            } else {
        //                Index blockRankZero = (BlockSize - blockRankOne) % BlockSize;
        //                return rankZero < blockRankZero;
        //            }
        //        };
        //        return std::upper_bound(blockRanks.ptr + lower, blockRanks.ptr + upper, rank, compareRank) - blockRanks.ptr;
        //        // TODO: ...until here
        // [[likely]] because for iid bits with probability p for '1', the expected value of upper - lower is 1 / p, ie 2 for 50% ones
        if (upper - lower <= linearFallbackSize) [[likely]] {
            // changing this from <= to < everywhere makes select 20% slower for no discernible reason, which is
            // especially weird because other choices of linearFallbackSize don't seem to matter much
            return selectBlockLinearly<IsOne>(lower, upper, rank);
        }
        // First, assume that ones are uniformly distributed. This may be completely wrong, in which case binary search
        // with that assumption could lead to linear running time, so quickly fall back to normal binary search.

        // Not really in danger of overflowing because rank < BlockSize and usually, BlockSize <= 256, upper < 2^48.
        Index expectedOffset = (upper - lower) * rank / BlockSize;
        ADS_ASSUME(expectedOffset >= 0);
        ADS_ASSUME(expectedOffset < upper - lower); // strict because rank < BlockSize, so use range [lower, upper] from now on
        Index mid = lower + expectedOffset;
        ADS_ASSUME(mid >= lower);
        ADS_ASSUME(mid < upper);
        Index midRank = getRank<IsOne>(mid);
        if (midRank <= rank) {
            lower = mid;       // don't add 1 in case mid + 1 == upper
            upper = upper - 1; // upper is now considered a valid part of the range
        } else {
            upper = mid;
        }
        ADS_ASSUME(lower <= upper);
        expectedOffset = (upper - lower) * rank / BlockSize;
        // get more pessimistic and try to get the expected position into a small interval
        // (instead of chopping away at only one end of the search interval)
        mid = lower + (expectedOffset * 3 + (upper - lower) / 2) / 4;
        ADS_ASSUME(mid >= lower);
        ADS_ASSUME(mid <= upper);
        midRank = getRank<IsOne>(mid);
        if (midRank <= rank) {
            lower = mid; // don't add 1 in case lower == upper
        } else {
            upper = mid;
        }
        while (upper - lower > linearFallbackSize) [[unlikely]] {
            mid = (lower + upper) / 2; // fallback to classical binary search
            ADS_ASSUME(mid >= lower);
            ADS_ASSUME(mid < upper);
            midRank = getRank<IsOne>(mid);
            if (midRank <= rank) {
                lower = mid + 1;
            } else {
                upper = mid;
            }
        }
        return selectBlockLinearly<IsOne>(lower, upper + 1, rank);
    }

    /// \brief The subrange [lower, upper) of block ranks must be in ascending order. Finds the first block in the range
    /// [lower, upper) where the rank is strictly larger than rank or returns upper if no such block exists.
    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectBlockLinearly(Index lower, Index upper, Index rank) const noexcept {
        ADS_ASSUME(rank >= 0);
        ADS_ASSUME(rank < BlockSize);
        ADS_ASSUME(lower >= 0);
        ADS_ASSUME(upper >= lower);
        ADS_ASSUME(upper - lower <= linearFallbackSize + 1);
        ADS_ASSUME(upper <= numBlocks());
        //        if constexpr (BlockSize <= 128) {
        //            const U64 mask = 0x8080'8080'8080'8080ull;
        //            U64 val = 42; // TODO: Load from U8* blockRanks + lower
        //            val |= mask;
        //            va -= rank * 0x1010'1010'1010'1010ull;
        //        }
        for (Index i = lower; i < upper; ++i) {
            if (getRank<IsOne>(i) > rank) { // TODO: Iterate over limbs instead of bytes?
                return i;
            }
        }
        return upper;
    }

    // TODO: Rename/remove
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankInBlock(const Limb* ADS_RESTRICT blockPtr) const noexcept {
        constexpr Index num256Blocks = BlockSize / 256; // encourage the compiler to unroll the "loop"
        Index res = 0;
        for (Index i = 0; i < num256Blocks; ++i) {
            res += u256Rank(blockPtr + 4 * i);
        }
        if constexpr (BlockSize % 256 != 0) {
            assert(BlockSize % 64 == 0);
            for (Index i = 0; i < (BlockSize % 256) / 64; ++i) {
                res += popcount(blockPtr[4 * num256Blocks + i]);
            }
        }
        return res;
    }
};


/// \brief Bitvector optimized for select queries. While TrivialBitvec may be faster, it also needs much more memory.
template<SupportedSelects SelectOps = SupportedSelects::BOTH, Index BlockSize = 256, typename BlockCountT = IntType<bytesNeededForIndexing(256)>>
using EfficientSelectBitvec
        = RecursiveBitvec<RecursiveBitvec<TrivialBitvec<U32, Operations::SELECT_ONLY>, SupportedSelects::ONE_ONLY>, SelectOps, BlockSize, BlockCountT>;


/// \brief The default bitvector type, which offers good performance across all operations. The same as EfficientSelectBitvec.
template<SupportedSelects SelectOps = SupportedSelects::BOTH, Index BlockSize = 256, typename BlockCountT = IntType<bytesNeededForIndexing(256)>>
using EfficientBitvec
        = RecursiveBitvec<RecursiveBitvec<TrivialBitvec<U32>, SupportedSelects::ONE_ONLY>, SelectOps, BlockSize, BlockCountT>;

static_assert(IsNormalBitvec<EfficientBitvec<>>);

} // namespace ads

#endif // ADS_RECURSIVE_BITVEC_HPP
