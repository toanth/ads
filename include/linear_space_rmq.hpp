#ifndef BITVECTOR_LINEAR_SPACE_RMQ_HPP
#define BITVECTOR_LINEAR_SPACE_RMQ_HPP

#include "common.hpp"
#include "nlogn_rmq.hpp"
#include <algorithm>
#include <cstdint>
namespace ads {

// TODO: For small sizes, use naive rmq instead of nlogn bitvecRmq?
template<typename T, Index BlockSize, typename InBlockIdx, typename BlockNumIdx, typename Comp>
struct NLogNBlockRmq : NLogNRmqOps<NLogNBlockRmq<T, BlockSize, InBlockIdx, BlockNumIdx, Comp>, Comp> {
    Index length = 0;
    const T* values = nullptr;
    InBlockIdx* minimumInBlock = nullptr;
    std::unique_ptr<BlockNumIdx[]> minima = nullptr; // TODO: Use shared allocation? Measure

    NLogNBlockRmq() = default;

    void build(const T* vals, Index len, InBlockIdx* blockMinima) noexcept {
        length = len;
        values = vals;
        minimumInBlock = blockMinima;
        minima = makeUniqueForOverwrite<BlockNumIdx>(minimaSize(length));
        this->init();
    }

    //    [[nodiscard]] T& getArrayElement(Index i) noexcept { return values[minimumInBlock[i]]; }
    [[nodiscard]] const T& getArrayElement(Index i) const noexcept { return values[i * BlockSize + minimumInBlock[i]]; }

    [[nodiscard]] Index getMinimum(Index i) const noexcept { return minima[i]; }
    void setMinimum(Index i, BlockNumIdx value) const noexcept { minima[i] = value; }
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
class LinearSpaceRMQ {
    static_assert(BlockSize % SubBlockSize == 0, "subblock size must divide block size");
    using InBlockIdx = IntType<bytesNeededForIndexing(BlockSize)>;
    using SubblockIndex = IntType<bytesNeededForIndexing(SubBlockStoresSubBlockIdx ? BlockSize / SubBlockSize : BlockSize)>;
    using BlockRmq = NLogNBlockRmq<T, BlockSize, InBlockIdx, BlockNumIdx, Comp>;
    // Stores the original array of values
    // size() elements; 64 * n bits for default params
    std::unique_ptr<T[]> values = nullptr;
    // O(n log n) space bitvecRmq structure over all blocks, which directly answers (sub)queries where lowre and upper
    // are multiples of the block size. size() / BlockSize elements (rounded up), each storing a logarithmic number of
    // indices; 32 * k log k (with k = n / 65536) bits by default; 64 * k log k, k = n / 256
    BlockRmq blockRmq = BlockRmq();
    // For each block, stores its mininum relative to the block start. Used by the bock bitvecRmq.
    // size() / BlockSize elements, rounded up; 16 * n / 65636 bits; 8 * n / 256 bits
    std::unique_ptr<InBlockIdx[]> minimumInBlock = nullptr;
    // For each sub-block, stores the exclusive prefix sum in its block.
    // size() /  SubBlockSize elements, possibly with up to (size() / BlockSize) additional elements due to rounding; 8 * n / 256 bits; 8 (maybe 4?) * n / 8 bits
    std::unique_ptr<SubblockIndex[]> subBlockPrefixMinima = nullptr;
    // For each sub-block, stores the exclusive suffix sum in its block.
    // size() /  SubBlockSize elements, possibly with up to (size() / BlockSize) additional elements due to rounding; 8 * n / 256 bits; 8 (maybe 4?) * n / 8 bits
    SubblockIndex* ADS_RESTRICT subBlockSuffixMinima = nullptr;
    // The number of entrues in the values array.
    Index length = 0;
    /// The total allocated size in bits.
    Index numAllocatedBits = 0;

    const T* ptrMin(const T* a, const T* b) const noexcept { return comp()(*a, *b) ? a : b; }

    Comp comp() const noexcept { return blockRmq.getComp(); }

public:
    using value_type = T;
    // TODO: Alternative idea: partition 256value blocks into 16value subblocks and store for each subblock the relative
    // index of its minimum
    //  in a 4 bit value, since there are 16 subblocks this requires 64 bits per block which may be stored directly
    //  after the block to reduce cache misses. These subblocks can speed up inner-block queries significantly. Also,
    //  subblock prefix and suffix arrays only need to store the subblock indices, which means each entry requires 4
    //  instead of 8 bit (it would probably make the most sense to increase their block size to 16 values as well)
    constexpr static Index blockSize = BlockSize;
    constexpr static Index subBlockSize = SubBlockSize;
    constexpr static const Index subBlocksPerBlock = roundUpDiv(BlockSize, SubBlockSize);
    constexpr static const char name[] = "linear space RMQ";

    LinearSpaceRMQ() = default;

    LinearSpaceRMQ(std::initializer_list<T> list) : LinearSpaceRMQ(Span<const T>(list.begin(), list.end())) {}

    explicit LinearSpaceRMQ(Span<const T> inputValues)
        : LinearSpaceRMQ(toUniquePtr(inputValues, roundUpDiv(inputValues.size(), SubBlockSize) * SubBlockSize),
                roundUpDiv(inputValues.size(), SubBlockSize) * SubBlockSize) {}

    explicit LinearSpaceRMQ(const std::vector<T>& inputValues) : LinearSpaceRMQ(Span<const T>(inputValues)) {}

    LinearSpaceRMQ(T* inputValues, Index length) : LinearSpaceRMQ(Span<const T>(inputValues, length)) {}

    /// \brief Take an unique_ptr and a size to avoid paying for the extra capacity() - size() space of a vector
    LinearSpaceRMQ(std::unique_ptr<T[]> inputValues, Index length)
        : values(std::move(inputValues)), blockRmq(), length(length) {
        assert(size() % SubBlockSize == 0); // The input array must have been extended by possibly values initialized elements
        const Index numBlocks = roundUpDiv(length, BlockSize);
        minimumInBlock = makeUniqueForOverwrite<InBlockIdx>(numBlocks);
        subBlockPrefixMinima = makeUniqueForOverwrite<SubblockIndex>(2 * subBlocksPerBlock * numBlocks);
        subBlockSuffixMinima = subBlockPrefixMinima.get() + subBlocksPerBlock * numBlocks;
        const T *prefixMinSoFar = nullptr, *suffixMinSoFar = nullptr;
        const T *beginOfBlock = nullptr, *endOfBlock = nullptr;
        const Index subBlockIdxQuotient = SubBlockStoresSubBlockIdx ? SubBlockSize : 1;

        for (Index i = 0; i < size(); i += SubBlockSize) {
            Index subBlockIdx = (i % BlockSize) / SubBlockSize;
            if (subBlockIdx == 0) {
                assert(!suffixMinSoFar || *suffixMinSoFar == *prefixMinSoFar);
                if (i > 0) {
                    assert(suffixMinSoFar >= beginOfBlock && suffixMinSoFar < endOfBlock);
                    minimumInBlock[i / BlockSize - 1] = suffixMinSoFar - beginOfBlock;
                }
                beginOfBlock = prefixMinSoFar = &values[i];
                endOfBlock = values.get() + std::min(i + BlockSize, size());
                suffixMinSoFar = endOfBlock - 1;
                subBlockPrefixMinima[i / SubBlockSize] = -1;
                //                Index blockLength = endOfBlock - beginOfBlock;
                subBlockSuffixMinima[(endOfBlock - values.get()) / SubBlockSize - 1] = -1; // blockLength / subBlockIdxQuotient - 1;
            }
            const T* iter = std::min_element(values.get() + i, values.get() + i + SubBlockSize);
            if (comp()(*iter, *prefixMinSoFar)) {
                prefixMinSoFar = iter;
            }
            if (subBlockIdx + 1 != subBlocksPerBlock) {
                subBlockPrefixMinima[i / SubBlockSize + 1] = (prefixMinSoFar - beginOfBlock) / subBlockIdxQuotient;
            }

            iter = std::min_element(endOfBlock - (subBlockIdx + 1) * SubBlockSize, endOfBlock - subBlockIdx * SubBlockSize);
            if (comp()(*iter, *suffixMinSoFar)) {
                suffixMinSoFar = iter;
            }
            if (subBlockIdx + 1 != subBlocksPerBlock) {
                subBlockSuffixMinima[(endOfBlock - values.get()) / SubBlockSize - subBlockIdx - 2]
                        = (suffixMinSoFar - beginOfBlock) / subBlockIdxQuotient;
            }
        }
        assert(size() == 0 || *suffixMinSoFar == *prefixMinSoFar);
        minimumInBlock[(size() - 1) / BlockSize] = suffixMinSoFar - beginOfBlock;
        assert(suffixMinSoFar >= beginOfBlock && suffixMinSoFar < endOfBlock);
        blockRmq.build(values.get(), numBlocks, minimumInBlock.get());
        numAllocatedBits = length * sizeof(T) + minimaSize(blockRmq.length) * sizeof(BlockNumIdx)
                           + numBlocks * sizeof(InBlockIdx) + 2 * subBlocksPerBlock * numBlocks * sizeof(SubblockIndex);
        numAllocatedBits *= 8;
    }


    [[nodiscard]] Index rmq(Index lower, Index upper) const {
        if (lower / BlockSize == (upper - 1) / BlockSize) {
            // TODO: Use information from minimumInBlock and possibly subblockPrefixMinimum / subblockSuffixMinimum
            return std::min_element(values.get() + lower, values.get() + upper) - values.get();
        }
        Index lowerCompleteBlock = roundUpDiv(lower, BlockSize);
        Index upperCompleteBlock = upper / BlockSize;
        assert(lowerCompleteBlock <= upperCompleteBlock);
        const T* minPtr = values.get() + lower;
        if (lowerCompleteBlock < upperCompleteBlock) {
            Index blockIdx = blockRmq(lowerCompleteBlock, upperCompleteBlock);
            minPtr = values.get() + blockIdx * BlockSize + minimumInBlock[blockIdx];
            assert(minPtr >= values.get() + lower && minPtr < values.get() + upper);
        }
        // TODO: Use minimumInBlock to sometimes speed up search in block prefix/suffix
        if (lower % BlockSize != 0) {
            const T* suffixMinPtr = values.get() + lower;
            Index lowerSubBlock = lower / SubBlockSize;
            if ((lowerSubBlock + 1) % subBlocksPerBlock != 0) {
                Index subBlockSuffixIdx = subBlockSuffixMinima[lowerSubBlock];
                const T* lowerBlockBeginPtr
                        = values.get() + (lowerCompleteBlock - 1) * BlockSize; // TODO: Reverse order in subBlockSuffixMinima?
                suffixMinPtr = lowerBlockBeginPtr + subBlockSuffixIdx;
                if constexpr (SubBlockStoresSubBlockIdx) {
                    suffixMinPtr = std::min_element(lowerBlockBeginPtr + subBlockSuffixIdx * SubBlockSize,
                            lowerBlockBeginPtr + (subBlockSuffixIdx + 1) * SubBlockSize);
                }
                assert(suffixMinPtr >= values.get() + lower && suffixMinPtr < values.get() + upper);
            }
            for (Index i = lower; i < (lowerSubBlock + 1) * SubBlockSize; ++i) {
                suffixMinPtr = ptrMin(values.get() + i, suffixMinPtr);
            }
            assert(suffixMinPtr >= values.get() + lower && suffixMinPtr < values.get() + upper);
            minPtr = ptrMin(minPtr, suffixMinPtr);
        }

        if (upper % BlockSize != 0) {
            Index upperSubBlock = upper / SubBlockSize;
            const T* prefixMinPtr = values.get() + upper - 1;
            if (upperSubBlock % subBlocksPerBlock != 0) {
                assert(subBlocksPerBlock < upperSubBlock && upperSubBlock <= size() / SubBlockSize);
                assert(upperSubBlock < subBlocksPerBlock * blockRmq.size());
                Index subBlockPrefixIdx = subBlockPrefixMinima[upperSubBlock];
                const T* upperBlockBeginPtr = values.get() + upperCompleteBlock * BlockSize;
                prefixMinPtr = upperBlockBeginPtr + subBlockPrefixIdx;

                if constexpr (SubBlockStoresSubBlockIdx) {
                    prefixMinPtr = std::min_element(upperBlockBeginPtr + subBlockPrefixIdx * SubBlockSize,
                            upperBlockBeginPtr + (subBlockPrefixIdx + 1) * SubBlockSize);
                }
                assert(prefixMinPtr >= values.get() + lower && prefixMinPtr < values.get() + upper);
            }
            for (Index i = upper - 1; i >= upperSubBlock * SubBlockSize; --i) {
                prefixMinPtr = ptrMin(values.get() + i, prefixMinPtr);
            }
            assert(prefixMinPtr >= values.get() + lower && prefixMinPtr < values.get() + upper);
            minPtr = ptrMin(minPtr, prefixMinPtr);
        }
        assert(minPtr >= values.get() + lower && minPtr < values.get() + upper);
        return minPtr - values.get();
    }

    Index operator()(Index lower, Index upper) const { return rmq(lower, upper); }

    [[nodiscard]] Index sizeInBits() const noexcept { return numAllocatedBits; }

    [[nodiscard]] Index size() const noexcept { return length; }
};

} // namespace ads

#endif // BITVECTOR_LINEAR_SPACE_RMQ_HPP
