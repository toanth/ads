#ifndef BITVECTOR_LINEAR_SPACE_RMQ_HPP
#define BITVECTOR_LINEAR_SPACE_RMQ_HPP

#include "common.hpp"
#include "nlogn_rmq.hpp"
#include <algorithm>
#include <cstdint>
namespace ads {

// TODO: For small sizes, use naive rmq instead of nlogn bitvecRmq?
template<typename T, Index BlockSize, typename InBlockIdx, typename BlockNumIdx, typename Comp>
struct [[nodiscard]] NLogNBlockRmq : NLogNRmqOps<NLogNBlockRmq<T, BlockSize, InBlockIdx, BlockNumIdx, Comp>, Comp> {
    Index length = 0;
    const T* values = nullptr;
    InBlockIdx* minimumInBlock = nullptr;
    View<BlockNumIdx> minima = View<BlockNumIdx>();
    //    std::unique_ptr<BlockNumIdx[]> minima = nullptr; // TODO: Use shared allocation? Measure
    constexpr NLogNBlockRmq() = default;

    ADS_CPP20_CONSTEXPR void build(const T* vals, Index len, InBlockIdx* blockMinima, Elem* minimaMemory) noexcept {
        length = len;
        values = vals;
        minimumInBlock = blockMinima;
        minima = View<BlockNumIdx>(minimaMemory, minimaSize(length));
        this->init();
    }

    //    [[nodiscard]] T& getArrayElement(Index i) noexcept { return values[minimumInBlock[i]]; }
    [[nodiscard]] ADS_CPP20_CONSTEXPR const T& getArrayElement(Index i) const noexcept {
        return values[i * BlockSize + minimumInBlock[i]];
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index getMinimum(Index i) const noexcept { return minima[i]; }
    ADS_CPP20_CONSTEXPR void setMinimum(Index i, BlockNumIdx value) const noexcept { minima.ptr[i] = value; }
};

// TODO: Don't copy values from an rvalue input vector into the RMQ struct, applies to all RMQ variants.
/// \brief Answers range minimum queries in "constant" time after linear preprocessing time while using linear space.
/// Additionally divides blocks into subblocks to answer in-block queries somewhat efficiently while not using too much
/// space. With the default template parameters, this class can only store arays of 2^48 elements, although that is easy
/// to change by setting BlockNumIdx to std::uint64_y. \tparam T \tparam BlockSize \tparam InBlockIdx \tparam
/// BlockNumIdx \tparam Comp
// template<typename T = Elem, Index BlockSize = 1 << 16, typename BlockNumIdx = std::uint32_t,
//          Index SubBlockSize = 1 << 8, typename Comp = std::less<>>
template<typename T = Elem, Index BlockSize = 1 << 8, typename BlockNumIdx = std::uint64_t, Index SubBlockSize = 1 << 3,
        typename Comp = std::less<>, bool SubBlockStoresSubBlockIdx = false>
class [[nodiscard]] LinearSpaceRMQ {
    static_assert(BlockSize % SubBlockSize == 0, "subblock size must divide block size");
    using InBlockIdx = IntType<bytesNeededForIndexing(BlockSize)>;
    using SubBlockIndex = IntType<bytesNeededForIndexing(SubBlockStoresSubBlockIdx ? BlockSize / SubBlockSize : BlockSize)>;
    using BlockRmq = NLogNBlockRmq<T, BlockSize, InBlockIdx, BlockNumIdx, Comp>;
    // The number of entries in the values array.
    Index length = 0;
    // Holds a pointer (and size) to all allocated memory; the destructor deallocates the memory. Doesn't call dtors.
    // Must come before any Views into the allocated data in the list of data members.
    Allocation<> allocation = Allocation();
    // Stores the original array of values
    // size() elements; 64 * n bits for default params
    View<T> values = View<T>();
    // O(n log n) space bitvecRmq structure over all blocks, which directly answers (sub)queries where lower and upper
    // are multiples of the block size. size() / BlockSize elements (rounded up), each storing a logarithmic number of
    // indices; 32 * k log k (with k = n / 65536) bits by default; 64 * k log k, k = n / 256
    BlockRmq blockRmq = BlockRmq();
    // For each block, stores its mininum relative to the block start. Used by the bock bitvecRmq.
    // size() / BlockSize elements, rounded up; 16 * n / 65636 bits; 8 * n / 256 bits
    View<InBlockIdx> minimumInBlock = View<InBlockIdx>();
    // For each sub-block, stores the exclusive prefix sum in its block.
    // size() /  SubBlockSize elements, possibly with up to (size() / BlockSize) additional elements due to rounding; 8 * n / 256 bits; 8 (maybe 4?) * n / 8 bits
    View<SubBlockIndex> subBlockPrefixMinima = View<SubBlockIndex>();
    // For each sub-block, stores the exclusive suffix sum in its block.
    // size() /  SubBlockSize elements, possibly with up to (size() / BlockSize) additional elements due to rounding; 8 * n / 256 bits; 8 (maybe 4?) * n / 8 bits
    View<SubBlockIndex> subBlockSuffixMinima = View<SubBlockIndex>();
    /// The total allocated size in bits.
    Index numAllocatedBits = 0;

    [[nodiscard]] constexpr const T* ptrMin(const T* a, const T* b) const noexcept { return comp()(*a, *b) ? a : b; }

    [[nodiscard]] constexpr Comp comp() const noexcept { return blockRmq.getComp(); }

public:
    using value_type = T;

    // TODO: Better idea: For each 256value block, use 128 bits to store a variant of the n log n RMQ: Partition
    // the block into 16 16value subblocks and build a structure similar to a mix of the rmmTree and the nlogn rmq over
    // subblocks: For each group of 2/4/8 consecutive subblocks, use 1 bit to store whether the left or right half
    // contains the smallest value. This uses 15 + 13 + 8 = 36 bits (or probably 16+16+8=40bits). For each subblock, use
    // 4 bits to store the index of its smallest element, using 4 * 16 = 64 bit. Additionally, it's possible to use 1
    // bit per subblock to sometimes speed up prefix/suffix queries within a subblock: If the smallest value in a
    // subblock is in the left half, store whether the smallest right-half value is less than the smallest value in the
    // left half after the smallest value in the left half; if the smallest subblock value is in the right half, store
    // whether the left half container a value smaller than all values in the right half before the smallest value
    // there. This often allows only looking at one half of a subblock, which corresponds to a cache line for 64 bit
    // values, and uses 16 bits of space. Also use 8 bits to store the index of the smallest value in the entire block,
    // using no more than 128 bits in total. In total, this can answer all rqm queries within the block (including
    // prefix and suffix queries) using 128 bits and reading at most 15 * 2 + 2 values, where 2 values are read to find
    // the smallest value within complete subblocks and typically, at most 8 consecutive and cacheline-aligned values
    // must be read to find the smallest value in a prefix/suffix. It should also be possible to navigate this metadata
    // without branching, using only bitwise operations. The metadata could be stored after the block values to improve
    // memory locality although a 256 value block already takes up exactly 1/2 page, so it's probably better to use a
    // separate array instead. A similar idea can be applied to superblocks of 256 blocks, except that it's better to
    // use more space to avoid having to read and compare block minima: Using the n log n rmq with base 2, block indices
    // take 8 bits, and there are 7 levels with approx. 256 entries each, totalling less than (probably exactly, due to
    // rounding up) 8 * 7 * 256 <= 2^14 bits, or 2^6 = 64 bits per block.
    constexpr static Index blockSize = BlockSize;
    constexpr static Index subBlockSize = SubBlockSize;
    constexpr static const Index subBlocksPerBlock = roundUpDiv(BlockSize, SubBlockSize);
    constexpr static const char name[] = "linear space RMQ";

    [[nodiscard]] static constexpr Index sizeInElems(Index length) noexcept {
        //        length = roundUpDiv(length, SubBlockSize) * SubBlockSize;
        assert(length % SubBlockSize == 0); // The input array must have been extended by possibly value initialized elements
        const Index numBlocks = roundUpDiv(length, BlockSize);
        Index sizeInBytes = length * sizeof(T) + minimaSize(numBlocks) * sizeof(BlockNumIdx)
                            + numBlocks * sizeof(InBlockIdx) + 2 * subBlocksPerBlock * numBlocks * sizeof(SubBlockIndex);
        return roundUpDiv(sizeInBytes, sizeof(Elem));
    }

    [[nodiscard]] static constexpr Index lengthInArray(Index length) noexcept {
        return roundUpDiv(length, BlockSize) * BlockSize;
    }

    constexpr LinearSpaceRMQ() = default;

    ADS_CPP20_CONSTEXPR LinearSpaceRMQ(std::initializer_list<T> list)
        : LinearSpaceRMQ(Span<const T>(list.begin(), list.end())) {}


    /// \brief Take an unique_ptr and a size to avoid paying for the extra capacity() - size() space of a vector
    ADS_CPP20_CONSTEXPR LinearSpaceRMQ(std::unique_ptr<T[]> inputValues, Index length)
        : LinearSpaceRMQ(Span<const T>(inputValues.get(), length)) {}

    explicit ADS_CPP20_CONSTEXPR LinearSpaceRMQ(const std::vector<T>& inputValues)
        : LinearSpaceRMQ(Span<const T>(inputValues)) {}

    ADS_CPP20_CONSTEXPR LinearSpaceRMQ(T* inputValues, Index length)
        : LinearSpaceRMQ(Span<const T>(inputValues, length)) {}

    explicit ADS_CPP20_CONSTEXPR LinearSpaceRMQ(Span<const T> inputValues)
        : length(lengthInArray(inputValues.size())), allocation(sizeInElems(length)), blockRmq() {
        assert(size() % SubBlockSize == 0); // The input array must have been extended, eg by value-initialized elements
        assert(allocation.size() == sizeInElems(length));
        numAllocatedBits = allocation.size() * sizeof(Elem) * 8;
        const Index numBlocks = roundUpDiv(length, BlockSize);
        values = View<T>(allocation.memory(), length);
        std::copy(inputValues.begin(), inputValues.end(), values.ptr);
        Index alreadyAllocated = roundUpDiv(length * sizeof(T), sizeof(Elem));
        Index newlyAllocated = roundUpDiv(numBlocks * sizeof(InBlockIdx), sizeof(Elem));
        minimumInBlock = View<SubBlockIndex>(allocation.memory() + alreadyAllocated, numBlocks);
        alreadyAllocated += newlyAllocated;
        newlyAllocated = roundUpDiv(subBlocksPerBlock * numBlocks * sizeof(SubBlockIndex), sizeof(Elem));
        subBlockPrefixMinima = View<SubBlockIndex>(allocation.memory() + alreadyAllocated, subBlocksPerBlock * numBlocks);
        alreadyAllocated += newlyAllocated;
        subBlockSuffixMinima = View<SubBlockIndex>(allocation.memory() + alreadyAllocated, subBlocksPerBlock * numBlocks);
        alreadyAllocated += newlyAllocated;
        const T *prefixMinSoFar = nullptr, *suffixMinSoFar = nullptr;
        const T *beginOfBlock = nullptr, *endOfBlock = nullptr;
        const Index subBlockIdxQuotient = SubBlockStoresSubBlockIdx ? SubBlockSize : 1;

        for (Index i = 0; i < size(); i += SubBlockSize) {
            Index subBlockIdx = (i % BlockSize) / SubBlockSize;
            assert(subBlockIdx >= 0 && subBlockIdx < subBlocksPerBlock);
            if (subBlockIdx == 0) {
                assert(!suffixMinSoFar || *suffixMinSoFar == *prefixMinSoFar);
                if (i > 0) {
                    assert(suffixMinSoFar >= beginOfBlock && suffixMinSoFar < endOfBlock);
                    minimumInBlock.ptr[i / BlockSize - 1] = suffixMinSoFar - beginOfBlock;
                }
                beginOfBlock = prefixMinSoFar = values.ptr + i;
                endOfBlock = values.ptr + i + BlockSize;
                assert(i + BlockSize <= size());
                suffixMinSoFar = endOfBlock - 1;
                subBlockPrefixMinima.ptr[i / SubBlockSize] = -1;
                subBlockSuffixMinima.ptr[(endOfBlock - values.ptr) / SubBlockSize - 1] = -1;
            }
            const T* iter = std::min_element(values.ptr + i, values.ptr + i + SubBlockSize);
            if (comp()(*iter, *prefixMinSoFar)) {
                prefixMinSoFar = iter;
            }
            if (subBlockIdx + 1 != subBlocksPerBlock) {
                subBlockPrefixMinima.ptr[i / SubBlockSize + 1] = (prefixMinSoFar - beginOfBlock) / subBlockIdxQuotient;
            }

            iter = std::min_element(endOfBlock - (subBlockIdx + 1) * SubBlockSize, endOfBlock - subBlockIdx * SubBlockSize);
            if (comp()(*iter, *suffixMinSoFar)) {
                suffixMinSoFar = iter;
            }
            if (subBlockIdx + 1 < subBlocksPerBlock) {
                subBlockSuffixMinima.ptr[(endOfBlock - values.ptr) / SubBlockSize - subBlockIdx - 2]
                        = (suffixMinSoFar - beginOfBlock) / subBlockIdxQuotient;
            }
        }
        assert(size() == 0 || *suffixMinSoFar == *prefixMinSoFar);
        minimumInBlock.ptr[(size() - 1) / BlockSize] = suffixMinSoFar - beginOfBlock;
        assert(suffixMinSoFar >= beginOfBlock && suffixMinSoFar < endOfBlock);
        blockRmq.build(values.ptr, numBlocks, minimumInBlock.ptr, allocation.memory() + alreadyAllocated);
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rmq(Index lower, Index upper) const {
        if (lower / BlockSize == (upper - 1) / BlockSize) {
            // TODO: Use information from minimumInBlock and possibly subblockPrefixMinimum / subblockSuffixMinimum
            return std::min_element(values.ptr + lower, values.ptr + upper) - values.ptr;
        }
        Index lowerCompleteBlock = roundUpDiv(lower, BlockSize);
        Index upperCompleteBlock = upper / BlockSize;
        assert(lowerCompleteBlock <= upperCompleteBlock);
        const T* minPtr = values.ptr + lower;
        if (lowerCompleteBlock < upperCompleteBlock) {
            Index blockIdx = blockRmq(lowerCompleteBlock, upperCompleteBlock);
            minPtr = values.ptr + blockIdx * BlockSize + minimumInBlock.ptr[blockIdx];
            assert(minPtr >= values.ptr + lower && minPtr < values.ptr + upper);
        }
        // TODO: Use minimumInBlock to sometimes speed up search in block prefix/suffix
        if (lower % BlockSize != 0) {
            const T* suffixMinPtr = values.ptr + lower;
            Index lowerSubBlock = lower / SubBlockSize;
            if ((lowerSubBlock + 1) % subBlocksPerBlock != 0) {
                Index subBlockSuffixIdx = subBlockSuffixMinima.ptr[lowerSubBlock];
                const T* lowerBlockBeginPtr
                        = values.ptr + (lowerCompleteBlock - 1) * BlockSize; // TODO: Reverse order in subBlockSuffixMinima?
                suffixMinPtr = lowerBlockBeginPtr + subBlockSuffixIdx;
                if constexpr (SubBlockStoresSubBlockIdx) {
                    suffixMinPtr = std::min_element(lowerBlockBeginPtr + subBlockSuffixIdx * SubBlockSize,
                            lowerBlockBeginPtr + (subBlockSuffixIdx + 1) * SubBlockSize);
                }
                assert(suffixMinPtr >= values.ptr + lower && suffixMinPtr < values.ptr + upper);
            }
            for (Index i = lower; i < (lowerSubBlock + 1) * SubBlockSize; ++i) {
                suffixMinPtr = ptrMin(values.ptr + i, suffixMinPtr);
            }
            assert(suffixMinPtr >= values.ptr + lower && suffixMinPtr < values.ptr + upper);
            minPtr = ptrMin(minPtr, suffixMinPtr);
        }

        if (upper % BlockSize != 0) {
            Index upperSubBlock = upper / SubBlockSize;
            const T* prefixMinPtr = values.ptr + upper - 1;
            if (upperSubBlock % subBlocksPerBlock != 0) {
                assert(subBlocksPerBlock < upperSubBlock && upperSubBlock <= size() / SubBlockSize);
                assert(upperSubBlock < subBlocksPerBlock * blockRmq.size());
                Index subBlockPrefixIdx = subBlockPrefixMinima[upperSubBlock];
                const T* upperBlockBeginPtr = values.ptr + upperCompleteBlock * BlockSize;
                prefixMinPtr = upperBlockBeginPtr + subBlockPrefixIdx;

                if constexpr (SubBlockStoresSubBlockIdx) {
                    prefixMinPtr = std::min_element(upperBlockBeginPtr + subBlockPrefixIdx * SubBlockSize,
                            upperBlockBeginPtr + (subBlockPrefixIdx + 1) * SubBlockSize);
                }
                assert(prefixMinPtr >= values.ptr + lower && prefixMinPtr < values.ptr + upper);
            }
            for (Index i = upper - 1; i >= upperSubBlock * SubBlockSize; --i) {
                prefixMinPtr = ptrMin(values.ptr + i, prefixMinPtr);
            }
            assert(prefixMinPtr >= values.ptr + lower && prefixMinPtr < values.ptr + upper);
            minPtr = ptrMin(minPtr, prefixMinPtr);
        }
        assert(minPtr >= values.ptr + lower && minPtr < values.ptr + upper);
        return minPtr - values.ptr;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index operator()(Index lower, Index upper) const { return rmq(lower, upper); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index sizeInBits() const noexcept { return numAllocatedBits; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return length; }
};

} // namespace ads

#endif // BITVECTOR_LINEAR_SPACE_RMQ_HPP
