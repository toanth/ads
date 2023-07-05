#ifndef ADS_RECURSIVE_BITVEC_HPP
#define ADS_RECURSIVE_BITVEC_HPP

#include "bitvec_base.hpp"
#include "trivial_bitvec.hpp"

namespace ads {

enum class SupportedSelectQueries { ZEROS, ONES, BOTH };

/// \brief Recursive Bitvector. Most of the later template arguments should not really be changed.
/// \tparam NestedBitvec The nested bitvector, of which two with approx. size() / BlockSize bits are stored.
/// \tparam Select Which select operations to support. Not supporting both saves one of two bitvectors storing approx. size()/BlockSize bits.
// TODO: Implement Space savings when Select is not BOTH?
/// \tparam BlockSize The number of bits in a block. This determines the BlockCountT type; smaller values use more space but can reduce query times.
/// \tparam BlockCountT The type used to represent block counts, must be able to hold values up to BlockSize - 1.
template<ADS_BITVEC_CONCEPT NestedBitvec, SupportedSelectQueries Select = SupportedSelectQueries::BOTH,
        Index BlockSize = 256, typename BlockCountT = IntType<bytesNeededForIndexing(BlockSize)>>
class RecursiveBitvec : public NormalBitvecBase<RecursiveBitvec<NestedBitvec, Select, BlockSize, BlockCountT>> {
    static_assert(sizeof(BlockCountT) >= bytesNeededForIndexing(BlockSize));
    static_assert(BlockSize >= 64); // TODO: Support smaller?

    using Base = NormalBitvecBase<RecursiveBitvec<NestedBitvec, Select, BlockSize, BlockCountT>>;
    friend Base;
    friend Base::Base;
    using Nested = NestedBitvec;

    /// \brief The nestedZeros Bitvector stores a bit for each block, which is one iff the block contains a zero at
    /// index i with rankZero(i) % BlockSize == 0. The nestedOnes Bitvector does the same for ones.
    Nested nestedZeros = Nested(); // store zeros bitvector before ones bitvector to keep the ones bitvector closer
    Nested nestedOnes = Nested();  // to the blockRanks array for better data locality on rank / selectOne
    /// \brief The blockRanks array stores the number of 1 at the end of the respective block, modulo BlockSize
    View<BlockCountT> blockRanks = View<BlockCountT>();
    /// \brief The actual bit sequence.
    View<Limb> vec = View<Limb>();
    /// \brief The logical number of bits in the bit sequence. The bitvector may internally allocate, set and access
    /// a greater number of bits that this value.
    Index numBits;

    constexpr static Index limbsInBlock = BlockSize / 64;

    constexpr static Index nestedSizeForBits(Index numBits) noexcept {
        // store a "1" after the actual nested bitvector to avoid a special case in select
        // also, the extra entry in the blockRanks array is used for numOnes(), in addition to special case avoidance
        return roundUpDiv(numBits, BlockSize) + 1;
    }

    ADS_CPP20_CONSTEXPR RecursiveBitvec(UninitializedTag, Index numBits, Limb* mem) noexcept
        : Base(numBits, mem),
          nestedZeros(Nested::uninitializedForSize(nestedSizeForBits(numBits), this->allocation.memory())),
          nestedOnes(Nested::uninitializedForSize(
                  nestedSizeForBits(numBits), this->allocation.memory() + nestedZeros.allocatedSizeInLimbs())),
          numBits(numBits) {
        numBits = roundUpDiv(numBits, BlockSize) * BlockSize; // don't deal with incomplete blocks
        Index nestedSize = nestedSizeForBits(numBits); // for the blockRanks array, the last entry can be used for numOnes()
        ADS_ASSUME(nestedSize == nestedOnes.size());
        mem = this->allocation.memory() + 2 * nestedZeros.allocatedSizeInLimbs();
        blockRanks = View<BlockCountT>(mem, nestedSize);
        vec = View<Limb>(this->allocation.memory() + this->allocation.size() - numBits / 64, numBits / 64);
        ADS_IF_CONSTEVAL {}
        else {
            ADS_ASSUME(std::greater_equal<void*>{}(vec.ptr, blockRanks.ptr + nestedSize));
        }
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
        Index nestedSizesInLimbs = Nested::allocatedSizeInLimbsForBits(numBlocks) * 2;
        Index arrSizeInLimbs = roundUpDiv(numBlocks * sizeof(BlockCountT), sizeof(Limb));
        Index bitSequenceSizeInLimbs = roundUpDiv(numBits, BlockSize) * limbsInBlock; // no incomplete blocks
        return nestedSizesInLimbs + arrSizeInLimbs + bitSequenceSizeInLimbs;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return numBits; }

    [[nodiscard]] constexpr Index numLimbs() const noexcept {
        ADS_ASSUME(vec.numT >= roundUpDiv(numBits, 64));
        return roundUpDiv(numBits, 64); // vec.numT rounds up towards full blocks
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBlocks() const noexcept { return nestedZeros.size() - 1; }

private:
    [[nodiscard]] ADS_CPP20_CONSTEXPR const Limb* getLimbArray() const noexcept { return vec.ptr; }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb* getLimbArray() noexcept { return vec.ptr; }

public:
    ADS_CPP20_CONSTEXPR void buildMetadata() noexcept {
        Index rankOne = 0;
        Index rankZero = 0;
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
            Index prevRankZero = rankZero;
            rankZero = (blockIdx + 1) * BlockSize - rankOne;
            nestedZeros.setBit(blockIdx, rankZero / BlockSize > prevRankZero / BlockSize);
        }
        ADS_ASSUME(rankOne <= size());
        blockRanks.ptr[numBlocks] = rankOne % BlockSize; // makes numOnes() simpler and avoids special cases for select(largeVal)
        nestedZeros.setBit(numBlocks);                   // more special case avoidance
        nestedOnes.setBit(numBlocks);
        nestedZeros.buildMetadata();
        nestedOnes.buildMetadata();
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


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectOneUnchecked(Index rank) const noexcept {
        ADS_ASSUME(rank >= 0);
        ADS_ASSUME(rank < numOnes());
        Index nestedSearchRank = rank / BlockSize;
        rank %= BlockSize;
        ADS_ASSUME(nestedSearchRank < nestedOnes.numOnes());
        Index nestedIdxEnd = nestedOnes.selectOne(nestedSearchRank);
        Index nestedIdx;
        if (nestedSearchRank == 0) [[unlikely]] {
            nestedIdx = 0; // if nestedIdx == nestedIdxEnd, std::upper_bound returns nestedIdx
        } else {
            nestedIdx = nestedOnes.selectOne(nestedSearchRank - 1);
            ADS_ASSUME(nestedIdxEnd > nestedIdx);
        }
        const BlockCountT* blockPtr = std::upper_bound(blockRanks.ptr + nestedIdx, blockRanks.ptr + nestedIdxEnd, rank);
        Index blockIdx = blockPtr - blockRanks.ptr;
        ADS_ASSUME(blockIdx >= nestedIdx);
        // blockRanks[blockIdx] contains the number of 1s relative to nestedSearchRank * BlockSize or (nestedSearchRank +
        // 1) * BlockSize, so we need to find the (blockRanks[blockIdx] - rank mod BlockSize)th bit counting from the right, 0 indexed
        Index rankInBlockFromRight = (blockRanks[blockIdx] + BlockSize - rank - 1) % BlockSize;
        ADS_ASSUME(rankInBlockFromRight >= 0);
        ADS_ASSUME(rankInBlockFromRight < rankInBlock(&this->getLimbRef(blockIdx * limbsInBlock)));
        return selectInBlockFromRight<true>(blockIdx, rankInBlockFromRight);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectZeroUnchecked(Index rankZero) const noexcept {
        ADS_ASSUME(rankZero >= 0);
        Index nestedSearchRank = rankZero / BlockSize;
        rankZero %= BlockSize;
        ADS_ASSUME(nestedSearchRank < nestedZeros.numOnes());
        Index nestedIdxEnd = nestedZeros.selectOne(nestedSearchRank);
        Index nestedIdx;
        if (nestedSearchRank == 0) [[unlikely]] {
            nestedIdx = 0;
        } else {
            nestedIdx = nestedZeros.selectOne(nestedSearchRank - 1); // a 1 represents a 256 * i th 0 with i > 0
            ADS_ASSUME(nestedIdxEnd > nestedIdx);
        }
        // The blockRank array stores rankOne() after the end of the block, modulo BlockSize.
        // rankZero((i + 1) * BlockSize) % BlockSize is ((i + 1) * BlockSize - blockRank[i]) % BlockSize, which is just
        // BlockSize - blockRank[i], except when the result would be BlockSize; it should be 0 instead.
        auto compareRank = [](Index rankZero, BlockCountT blockRankOne) {
            Index blockRankZero = (BlockSize - blockRankOne) % BlockSize;
            return rankZero < blockRankZero;
        };
        const BlockCountT* blockPtr
                = std::upper_bound(blockRanks.ptr + nestedIdx, blockRanks.ptr + nestedIdxEnd, rankZero, compareRank);
        Index blockIdx = blockPtr - blockRanks.ptr;
        ADS_ASSUME(blockIdx >= nestedIdx);
        ADS_ASSUME(blockIdx < numBlocks());
        // blockRanks[blockIdx] contains the number of 1s relative to nestedSearchRan * BlockSize or (nestedSearchRank +
        // 1) * BlockSize, so we need to find the (blockRans[blockIdx] - rankZero mod BlockSize)th bit counting from
        // the right, 0 indexed
        Index rankInBlockFromRight = (2 * BlockSize - blockRanks[blockIdx] - rankZero - 1) % BlockSize;
        if (blockIdx == numBlocks() - 1) [[unlikely]] {
            // TODO: Store complete blocks, make sure all limbs after numLimbs() are 0
            Index numIgnoredZeros = ((limbsInBlock - numLimbs() % limbsInBlock) % limbsInBlock) * 64;
            rankInBlockFromRight -= numIgnoredZeros;
        }
        ADS_ASSUME(rankInBlockFromRight >= 0); // TODO: Does this hold? explain why
        ADS_ASSUME(rankInBlockFromRight < BlockSize);
        return selectInBlockFromRight<false>(blockIdx, rankInBlockFromRight);
    }

private:
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
template<SupportedSelectQueries Select = SupportedSelectQueries::BOTH, Index BlockSize = 256, typename BlockCountT = IntType<bytesNeededForIndexing(256)>>
using EfficientSelectBitvec
        = RecursiveBitvec<RecursiveBitvec<TrivialBitvec<U32, Operations::SELECT_ONLY>, SupportedSelectQueries::ONES>, Select, BlockSize, BlockCountT>;

/// \brief The default bitvector type, which offers good performance across all operations. The same as EfficientSelectBitvec.
template<SupportedSelectQueries Select = SupportedSelectQueries::BOTH, Index BlockSize = 256, typename BlockCountT = IntType<bytesNeededForIndexing(256)>>
using EfficientBitvec = EfficientSelectBitvec<Select, BlockSize, BlockCountT>;

static_assert(IsNormalBitvec<EfficientBitvec<>>);

} // namespace ads

#endif // ADS_RECURSIVE_BITVEC_HPP
