#ifndef ADS_RECURSIVE_BITVEC_HPP
#define ADS_RECURSIVE_BITVEC_HPP

#include "base/normal_bitvec.hpp"
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
class RecursiveBitvec
    : public NormalBitvecBase<RecursiveBitvec<NestedBitvec, SelectOps, BlockSize, BlockCountT>, BlockSize / 64> {
    static_assert(sizeof(BlockCountT) >= bytesNeededForIndexing(BlockSize));
    static_assert(BlockSize >= 64);
    static_assert(BlockSize % 64 == 0);

    using Base = NormalBitvecBase<RecursiveBitvec<NestedBitvec, SelectOps, BlockSize, BlockCountT>, BlockSize / 64>;
    friend Base;
    friend Base::Base;
    using Nested = NestedBitvec;

    /// \brief The nestedZeros Bitvector stores a bit for each block b, which is one iff the block contains a zero at
    /// index i in [b * BlockSize, (b+1) * BlockSize) with rankZero(i) % BlockSize == 0 and rankZero(i) > 0.
    /// The nestedOnes Bitvector does the same for ones.
    Nested nestedZeros = Nested(); // store zeros bitvector before ones bitvector to keep the ones bitvector closer
    Nested nestedOnes = Nested();  // to the blockRanks array for better data locality on rank / selectOne
    /// \brief The blockRanks array stores the number of 1 at the end of the respective block, modulo BlockSize.
    Array<BlockCountT> blockRanks = Array<BlockCountT>();
    /// \brief The actual bit sequence.
    Array<CacheLine> vec = Array<CacheLine>();

    constexpr static Index limbsInBlock = BlockSize / 64;

    constexpr static Index nestedSizeForBits(Index numBits) noexcept {
        // append a "1" to the nested bitvector to avoid a special case in select
        // also, the extra entry in the blockRanks array is additionally used for numOnes()
        return roundUpDiv(numBits, BlockSize) + 1;
    }

    constexpr static Index bitSequenceSizeInCacheLines(Index numBits) noexcept {
        // round up to full blocks and full cache lines
        return roundUpDiv(numBits, CACHELINE_SIZE_BYTES * 8);
    }

    ADS_CPP20_CONSTEXPR RecursiveBitvec(UninitializedTag, Index numBits, CacheLine* mem) : Base(numBits, mem) {
        mem = this->allocation.memory();
        Index numCacheLinesInBitSequence = bitSequenceSizeInCacheLines(numBits);
        vec = Array<CacheLine>(mem, numCacheLinesInBitSequence);
        ADS_ASSUME_ALIGNED(vec.ptr, CACHELINE_SIZE_BYTES);
        mem += numCacheLinesInBitSequence;
        Index nestedSize = nestedSizeForBits(numBits); // for the blockRanks array, the last entry can be used for numOnes()
        if constexpr (SelectOps != SupportedSelects::ONE_ONLY) {
            nestedZeros = Nested::uninitializedForSize(nestedSize, mem);
            mem += nestedZeros.template allocatedSizeIn<Group::CacheLine>();
            ADS_ASSUME(nestedSize == nestedZeros.size());
        }
        nestedOnes = Nested::uninitializedForSize(nestedSize, mem);
        mem += nestedOnes.template allocatedSizeIn<Group::CacheLine>(); // TODO: For the TrivialBitvec, it's unnecessary to waste the remaining cache line
        blockRanks = Array<BlockCountT>(mem, nestedSize);
        ADS_ASSUME(nestedSize == nestedOnes.size());
        ADS_ASSUME(mem + roundUpDiv(blockRanks.sizeInBytes(), CACHELINE_SIZE_BYTES)
                   == this->allocation.memory() + this->allocation.sizeInTs());
    }


public:
    constexpr RecursiveBitvec() noexcept = default;

    explicit constexpr RecursiveBitvec(Index numBits, Limb fill, CacheLine* mem = nullptr)
        : RecursiveBitvec(UninitializedTag{}, numBits, mem) {
        this->fill(fill);
    }

    explicit ADS_CPP20_CONSTEXPR RecursiveBitvec(Span<const Limb> limbs, CacheLine* mem = nullptr) noexcept
        : RecursiveBitvec(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR RecursiveBitvec(Span<const Limb> limbs, Index numBits, CacheLine* mem = nullptr) noexcept
        : RecursiveBitvec(UninitializedTag{}, numBits, mem) {
        this->copyFrom(limbs);
    }

    explicit ADS_CPP20_CONSTEXPR RecursiveBitvec(std::string_view str, Index base = 2, CacheLine* mem = nullptr) noexcept
        : RecursiveBitvec(UninitializedTag{}, str.size() * intLog2(base), mem) {
        this->initFromStr(str, base);
    }

    [[nodiscard]] static constexpr Index allocatedSizeInBytesForBits(Index numBits) noexcept {
        Index numBlocks = nestedSizeForBits(numBits);
        const Index factor = (SelectOps == SupportedSelects::BOTH ? 2 : 1);
        Index nestedSizesInBytes = roundUpTo(Nested::allocatedSizeInBytesForBits(numBlocks), CACHELINE_SIZE_BYTES) * factor;
        Index arrSizeInBytes = roundUpTo(numBlocks * sizeof(BlockCountT), CACHELINE_SIZE_BYTES);
        static_assert(RecursiveBitvec::requiredAlignment() % Nested::requiredAlignment() == 0);
        static_assert(alignof(Limb) % alignof(BlockCountT) == 0);
        return bitSequenceSizeInCacheLines(numBits) * CACHELINE_SIZE_BYTES + nestedSizesInBytes + arrSizeInBytes;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBlocks() const noexcept { return nestedOnes.size() - 1; }


    [[nodiscard]] static constexpr Index numLimbsInBlock() noexcept { return BlockSize / 64; }

    [[nodiscard]] static constexpr Index blockSize() noexcept { return BlockSize; }



private:
    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<const CacheLine> getCacheLineArray() const noexcept {
        return {vec.ptr, vec.numT};
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<CacheLine> getCacheLineArray() noexcept { return {vec.ptr, vec.numT}; }

    //    ADS_CPP20_CONSTEXPR void completeWithZeros() noexcept { /*do nothing*/
    //    }

public:
    ADS_CPP20_CONSTEXPR void buildMetadata() noexcept {
        [[maybe_unused]] Index rankOne = 0;
        [[maybe_unused]] Index rankZero = 0;
        Index numBlocks = this->numBlocks();
        assert(numBlocks == this->template sizeIn<Group::Block>());
        this->completeWithZeros();

        for (Index blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
            Index prevRankOne = rankOne;
            rankOne += rankInBlock(blockIdx);
            nestedOnes.setBit(blockIdx, rankOne / BlockSize > prevRankOne / BlockSize);
            blockRanks.ptr[blockIdx] = rankOne % BlockSize;
            if constexpr (SelectOps != SupportedSelects::ONE_ONLY) {
                Index prevRankZero = rankZero;
                rankZero = (blockIdx + 1) * BlockSize - rankOne;
                nestedZeros.setBit(blockIdx, rankZero / BlockSize > prevRankZero / BlockSize);
            }
        }
        ADS_ASSUME(rankOne <= this->size());
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
        Index blockEndLimbIdx = (blockIdx + 1) * numLimbsInBlock();
        // count the number of 1s between idx and the end of the block. This doesn't require looking at
        // the previous block rank, which helps a lot when blockIdx is 0.
        for (Index i = blockEndLimbIdx - 1; i > limbIdx; --i) { // TODO: Unroll? SIMD? Implement same for rankZero?
            res -= popcount(this->getLimb(i));
        }
        const Limb mask = (Limb(-1) << (idx % 64));
        return res - popcount(this->getLimb(limbIdx) & mask);
    }

    // Don't iterate over the blockRanks twice and don't call nested selectOne twice unless necessary.
    // Cf. rankOneUnchecked() below. This function exists because it fits the pattern in which Elias-Fano calls
    // selectOne and the pattern in which this bitvector calls the nested bitvector's selectOne function, so for a RecursiveBitvec<RecursiveBitvec<TrivialBitvec>>>,
    // this turns 14 selectOne()s per predecessor into usually no more than 3 selectOneAndPrev calls overall.
    [[nodiscard]] [[using gnu: hot, pure]] ADS_CPP20_CONSTEXPR std::pair<Index, Index> selectOneAndPrevOne(
            Index rankOfSecond) const noexcept {
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
        ADS_ASSUME(rankInBlockFromRight < rankInBlock(blockIdx));
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
        ADS_ASSUME(rankInBlockFromRight < rankInBlock(blockIdx));
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
        //        if (blockIdx == numBlocks() - 1) [[unlikely]] {
        //            // TODO: Store complete blocks, make sure all limbs after numLimbs() are 0
        //            Index numIgnoredZeros = ((limbsInBlock - numLimbs() % limbsInBlock) % limbsInBlock) * 64;
        //            rankInBlockFromRight -= numIgnoredZeros;
        //        }
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
        Index blockEnd = (blockIdx + 1) * limbsInBlock;
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
        ADS_ASSUME(rankOfRight + 1 < BlockSize);        // the left value's rank must also be in the same block
        Index blockStart = blockIdx * limbsInBlock;
        Index blockEnd = (blockIdx + 1) * limbsInBlock; // this bitvector always stores full blocks
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

    /// \brief The subrange [lower, upper) of block ranks must be in ascending order. Finds the first block in the
    /// range [lower, upper) where the rank is strictly larger than rank or returns upper if no such block exists.
    template<bool IsOne>
    [[nodiscard]] [[using gnu: hot, pure]] ADS_CPP20_CONSTEXPR Index selectBlock(Index lower, Index upper, Index rank) const noexcept {
        ADS_ASSUME(rank >= 0);
        ADS_ASSUME(rank < BlockSize);
        ADS_ASSUME(lower >= 0);
        ADS_ASSUME(upper >= lower);
        ADS_ASSUME(upper <= numBlocks());
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

    /// \brief The subrange [lower, upper) of block ranks must be in ascending order. Finds the first block in the
    /// range [lower, upper) where the rank is strictly larger than rank or returns upper if no such block exists.
    template<bool IsOne>
    [[nodiscard]] [[using gnu: hot, pure]] ADS_CPP20_CONSTEXPR Index selectBlockLinearly(
            Index lower, Index upper, Index rank) const noexcept {
        ADS_ASSUME(rank >= 0);
        ADS_ASSUME(rank < BlockSize);
        ADS_ASSUME(lower >= 0);
        ADS_ASSUME(upper >= lower);
        ADS_ASSUME(upper - lower <= linearFallbackSize + 1);
        ADS_ASSUME(upper <= this->numAccessibleBlocks());
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

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankInBlock(Index blockIdx) const noexcept {
        Index res = 0;
        static_assert(BlockSize % 64 == 0);
        if constexpr (BlockSize % 256 == 0) {
            for (Index i = 0; i < BlockSize / 256; ++i) {
                res += u256Rank(&this->getLimbRef(blockIdx * this->numLimbsInBlock() + 4 * i));
            }
        } else {
            for (Index i = 0; i < BlockSize / 64; ++i) {
                Limb limb = this->getLimb(blockIdx * this->numLimbsInBlock() + i);
                res += popcount(limb);
            }
        }
        return res;
    }
};


/// \brief Bitvector optimized for select queries. While TrivialBitvec may be faster, it also needs much more memory.
template<SupportedSelects SelectOps = SupportedSelects::BOTH, Index BlockSize = 256, typename BlockCountT = IntType<bytesNeededForIndexing(256)>>
using EfficientSelectBitvec
        = RecursiveBitvec<RecursiveBitvec<TrivialBitvec<U32, Operations::SELECT_ONLY>, SupportedSelects::ONE_ONLY, 64>, SelectOps, BlockSize, BlockCountT>;


/// \brief The default bitvector type, which offers good performance across all operations. The same as EfficientSelectBitvec.
template<SupportedSelects SelectOps = SupportedSelects::BOTH, Index BlockSize = 256, typename BlockCountT = IntType<bytesNeededForIndexing(256)>>
using EfficientBitvec
        = RecursiveBitvec<RecursiveBitvec<TrivialBitvec<U32>, SupportedSelects::ONE_ONLY, 64>, SelectOps, BlockSize, BlockCountT>;

static_assert(IsNormalBitvec<EfficientBitvec<>>);

} // namespace ads

#endif // ADS_RECURSIVE_BITVEC_HPP
