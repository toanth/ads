
#ifndef ADS_SUCCINCT_RMQ_HPP
#define ADS_SUCCINCT_RMQ_HPP

#include "bitvector/classical_rank_bitvec.hpp"
#include "bitvector/recursive_bitvec.hpp"
#include "common.hpp"
#include <algorithm>
#include <stack>
namespace ads {

/// This file implements the succinct rmq data structure from https://algo2.iti.kit.edu/download/rmq.pdf,
/// using the range min-max tree from https://arxiv.org/pdf/1601.06939.pdf for solving the excess rmq problem on bitvectors.

// TODO: Flip openParen and closeParen, use SupportedSelects::ONE_ONLY for the bitvector
constexpr static bool openParen = true;
constexpr static bool closeParen = false;

using DefaultBitvec = EfficientBitvec<>;

template<typename T = Limb, Index BlockSize = 512, ADS_NORMAL_BITVEC_CONCEPT Bitvec = DefaultBitvec>
struct [[nodiscard]] RangeMinMaxTree {
    static_assert(BlockSize * 2 <= std::numeric_limits<T>::max());
    static_assert(BlockSize % 64 == 0);

    using value_type = T;

    static constexpr Index blockSize = BlockSize;

    Index numLeaves = 0;
    Index size = 0;
    Array<T> rmmArr = Array<T>();
    Bitvec bv;

    constexpr RangeMinMaxTree() = default;

    ADS_CPP20_CONSTEXPR ~RangeMinMaxTree() noexcept = default;

    ADS_CPP20_CONSTEXPR RangeMinMaxTree(RangeMinMaxTree&&) noexcept = default;

    ADS_CPP20_CONSTEXPR RangeMinMaxTree& operator=(RangeMinMaxTree&&) noexcept = default;
    ADS_CPP20_CONSTEXPR RangeMinMaxTree& operator=(const RangeMinMaxTree&) noexcept = delete;

    [[nodiscard]] constexpr static Index numLeavesForBits(Index numBits) noexcept {
        return roundUpDiv(numBits, BlockSize);
    }

    [[nodiscard]] constexpr static Index numNodesInArray(Index numLeaves) noexcept {
        return numLeaves + (Limb(1) << roundUpLog2(Limb(numLeaves))); // the missing -1 makes this a 1-indexed heap
    }

    [[nodiscard]] constexpr Index leafIdxInArr(Index leafNum) const noexcept { return size - numLeaves + leafNum; }

    /// The reverse of leafIdxInArr
    [[nodiscard]] constexpr Index leafNum(Index leafIdxInArr) const noexcept {
        assert(leafIdxInArr >= size - numLeaves && leafIdxInArr < size);
        return leafIdxInArr - (size - numLeaves);
    }

    [[nodiscard]] constexpr bool isRightChild(Index node) const noexcept { return node % 2 == 1; }

    // TODO: Different implementation: Store for each leaf the position of its minimum and the minimum excess in block
    // (>= -2^BlockLength, <= 2^BlockLength), which can be done with 16 bits per leaf if the block size is no greater
    // than 256 (min. excess in block can be relative to next block to keep numbers small). Implement and measure if
    // faster (for the same block size, this additionally uses 16 bits per leaf, which is roughly 8 bits per node but
    // may cause cache misses)

    template<typename U>
    explicit ADS_CPP20_CONSTEXPR RangeMinMaxTree(Bitvec&& bitvector, U* rmmArrPtr) noexcept
        : numLeaves(numLeavesForBits(bitvector.size())), size(numNodesInArray(numLeaves)), rmmArr(rmmArrPtr, size),
          bv(std::move(bitvector)) {
        numLeaves = numLeavesForBits(bv.size());
        assert(numLeaves > 0);
        size = numNodesInArray(numLeaves);
        assert(leafIdxInArr(0) % 2 == 0 || leafIdxInArr(0) == 1);
        Index minExcessInBlock = 1;
        Index excess = 0;
        for (Index i = 0; i < bv.size(); ++i) {
            if (i > 0 && i % BlockSize == 0) {
                Index l = leafIdxInArr(i / BlockSize - 1);
                rmmArr.ptr[l] = minExcessInBlock;
                while (isRightChild(l)) {
                    rmmArr.ptr[l / 2] = std::min(rmmArr[l - 1], rmmArr[l]);
                    l /= 2;
                }
                minExcessInBlock = excess + 1;
            }
            bool bit = bv.getBit(i);
            if (bit == openParen) { // TODO: There's probably a better way to do this than to read each bit individually - read bytes and use lookup table?
                ++excess;
            } else {
                --excess;
                minExcessInBlock = std::min(minExcessInBlock, excess);
            }
        }
        assert(minExcessInBlock == 0);
        Index v = leafIdxInArr(numLeaves - 1);
        rmmArr.ptr[v] = 0;
        for (; v > 0; v /= 2) {
            if (isRightChild(v)) {
                if (v / 2 == 0) { // TODO: Remove
                    rmmArr.ptr[v / 2] = rmmArr[v - 1];
                }
                rmmArr.ptr[v / 2] = std::min(rmmArr[v], rmmArr[v - 1]);
            } else {
                rmmArr.ptr[v / 2] = rmmArr[v];
            }
        }
    }
    void printTree() const {
        std::cout << " tree stored in array with " << size << " values, " << numLeaves << " leaves:" << std::endl;
        for (Index i = 1; i < this->size; ++i) {
            if (i == this->leafIdxInArr(0)) {
                std::cout << "(leafs) ";
            }
            std::cout << this->rmmArr[i] << " ";
            if (this->isRightChild(i)) {
                std::cout << ' ';
            }
            if ((i & (i + 1)) == 0) {
                std::cout << std::endl;
            }
        }
        std::cout << "\tend of tree" << std::endl;
    }

    struct MinRes {
        Index pos;
        Index minExcess;
    };

    static constexpr MinRes noRes = {std::numeric_limits<Index>::max(), std::numeric_limits<Index>::max()};

    [[nodiscard]] ADS_CPP20_CONSTEXPR MinRes findMinInBlockImpl(Index lower, Index upper, Index& excessSoFar) const noexcept {
        if (lower == upper) {
            return noRes;
        }
        Index excess = excessSoFar;
        Index minExcess = excess + 1;
        Index minPos = lower;
        for (Index i = lower; i < upper; ++i) {
            if (bv.getBit(i) == openParen) {
                ++excess;
            } else if (--excess < minExcess) { // choose leftmost index
                minExcess = excess;
                minPos = i;
            }
        }
        excessSoFar = excess;
        return {minPos, minExcess};
    }

    // Idea: find min extent within block and use that together with the heap's stored value of the global min extent in this block
    // to figure out the extent at the start of the block, which allows calculating the global min extent within the query range
    [[nodiscard]] ADS_CPP20_CONSTEXPR MinRes findMinInBlock(Index lower, Index upper) const noexcept {
        Index blockBegin = lower / BlockSize * BlockSize;
        Index blockEnd = blockBegin + BlockSize;
        blockEnd = std::min(blockEnd, bv.size());
        upper = std::min(upper, bv.size());
        assert(lower <= upper);
        if (upper == lower) {
            return noRes;
        }
        Index inBlockExcess = 0;
        Index minExcessInBlock = findMinInBlockImpl(blockBegin, lower, inBlockExcess).minExcess;
        MinRes res = findMinInBlockImpl(lower, upper, inBlockExcess);
        assert(res.pos >= lower && res.pos < upper);
        minExcessInBlock = std::min(minExcessInBlock, res.minExcess);
        Index minExcessAfterQuery = findMinInBlockImpl(upper, blockEnd, inBlockExcess).minExcess;
        minExcessInBlock = std::min(minExcessInBlock, minExcessAfterQuery);
        Index globalOffset = rmmArr[leafIdxInArr(lower / BlockSize)] - minExcessInBlock;
        return {res.pos, res.minExcess + globalOffset};
    }

    /// \brief Find the minimum element in the half-open range of blocks [\p lowerBlockIdx, \p upperBlockIdx)
    /// \param lowerBlockIdx
    /// \param upperBlockIdx
    /// \return
    [[nodiscard]] ADS_CPP20_CONSTEXPR MinRes findMinInTree(Index lowerBlockIdx, Index upperBlockIdx) const noexcept {
        assert(lowerBlockIdx <= upperBlockIdx + 1);
        if (upperBlockIdx <= lowerBlockIdx) {
            return noRes;
        } else if (lowerBlockIdx + 1 == upperBlockIdx) {
            return findMinInBlock(lowerBlockIdx * BlockSize, upperBlockIdx * BlockSize);
        }
        Index lower = leafIdxInArr(lowerBlockIdx);
        Index upper = leafIdxInArr(upperBlockIdx - 1);
        assert(intLog2(lower) == intLog2(upper));
        Index h = intLog2(lower ^ upper) + 1;
        Index lca = lower >> h;
        Index v = lower;
        // TODO: This can probably be done with some bit fiddling instead of a loop (or look at rmmArr and don't set leftMinPos unless necessary), also for rightMinPos
        while (!isRightChild(v) && v / 2 > lca) {
            assert(v > lca);
            v /= 2;
        }
        Index leftMinPos = v;
        //        v /= 2;
        while (v / 2 > lca) {
            if (!isRightChild(v)) {
                if (rmmArr[v + 1] < rmmArr[leftMinPos]) {
                    leftMinPos = v + 1;
                }
            }
            v /= 2;
        }
        v = upper;
        while (isRightChild(v) && v / 2 > lca) {
            assert(v > lca);
            v /= 2;
        }
        Index rightMinPos = v;
        //        v /= 2;
        while (v / 2 > lca) {
            if (isRightChild(v)) {
                if (rmmArr[v - 1] <= rmmArr[rightMinPos]) {
                    rightMinPos = v - 1;
                }
            }
            v /= 2;
        }
        Index minPos = rmmArr[rightMinPos] < rmmArr[leftMinPos] ? rightMinPos : leftMinPos; // choose left on ties
        while (minPos < lower) {
            minPos *= 2;
            assert(minPos <= upper);
            if (rmmArr[minPos + 1] < rmmArr[minPos]) {
                ++minPos;
            }
            assert(minPos <= upper);
        }
        Index i = leafNum(minPos);
        assert(std::min_element(rmmArr.ptr + lower, rmmArr.ptr + upper + 1) - rmmArr.ptr == minPos);
        return findMinInBlock(i * BlockSize, (i + 1) * BlockSize);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR MinRes rmqImpl(Index lower, Index upper) const noexcept {
        assert(lower < upper);
        Index nextLowerBlockIdx = roundUpDiv(lower, BlockSize);
        Index prevUpperBlockIdx = upper / BlockSize;
        if (prevUpperBlockIdx + 1 == nextLowerBlockIdx) {
            return findMinInBlock(lower, upper);
        }
        MinRes lowerCandidate = findMinInBlock(lower, nextLowerBlockIdx * BlockSize);
        assert((lowerCandidate.pos >= lower && lowerCandidate.pos < upper) || lowerCandidate.pos == noRes.pos);
        MinRes treeCandidate = findMinInTree(nextLowerBlockIdx, prevUpperBlockIdx);
        assert((treeCandidate.pos >= lower && treeCandidate.pos < upper) || treeCandidate.pos == noRes.pos);
        MinRes lowerOrTreeCandidate = treeCandidate.minExcess < lowerCandidate.minExcess ? treeCandidate : lowerCandidate; // take lower when tied
        MinRes upperCandidate = findMinInBlock(prevUpperBlockIdx * BlockSize, upper);
        assert((upperCandidate.pos >= lower && upperCandidate.pos < upper) || upperCandidate.pos == noRes.pos);
        return upperCandidate.minExcess < lowerOrTreeCandidate.minExcess ? upperCandidate : lowerOrTreeCandidate; // don't take upper when tied
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index bitvecRmq(Index lower, Index upper) const noexcept {
        return rmqImpl(lower, upper).pos;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const Bitvec& getBitvector() const noexcept { return bv; }
};


template<typename Bitvec = DefaultBitvec>
class [[nodiscard]] SuccinctRMQ {

    using RmmTree = RangeMinMaxTree<>;

    Allocation<> allocation = Allocation();
    RmmTree rmmTree = RmmTree();
    Index length = 0;

    template<typename Container>
    [[nodiscard]] static ADS_CPP20_CONSTEXPR Index numAllocatedBytes(const Container& container) noexcept {
        Index numBits = container.size() * 2 + 2;
        return Bitvec::template allocatedSizeForBitsIn<Group::CacheLine>(numBits) * CACHELINE_SIZE_BYTES
               + roundUpTo(sizeof(RmmTree::value_type) * RmmTree::numNodesInArray(RmmTree::numLeavesForBits(numBits)),
                       CACHELINE_SIZE_BYTES);
    }

    template<typename Container, typename Comp>
    ADS_CPP20_CONSTEXPR static Bitvec createDfuds(const Container& container, Comp comp, CacheLine* mem) noexcept {
        Bitvec dfuds = Bitvec::uninitializedForSize(container.size() * 2 + 2, mem);
        Index bvIdx = dfuds.size();
        std::vector<Index> stack;
        for (Index i = container.size() - 1; i + 1 > 0; --i) {
            assert(bvIdx > 1);
            dfuds.setBit(--bvIdx, closeParen);
            while (!stack.empty() && comp(container[i], container[stack.back()])) {
                assert(bvIdx > 1);
                dfuds.setBit(--bvIdx, openParen);
                stack.pop_back();
            }
            stack.push_back(i);
        }
        dfuds.setBit(--bvIdx, closeParen);
        while (!stack.empty()) {
            stack.pop_back();
            dfuds.setBit(--bvIdx, openParen);
        }
        assert(bvIdx == 1);
        dfuds.setBit(0, openParen);
        dfuds.buildMetadata();
        return dfuds;
    }

    template<typename Container, typename Comp>
    ADS_CPP20_CONSTEXPR RmmTree createRmmTree(const Container& container, Comp comp) const noexcept {
        CacheLine* ptr = allocation.memory();
        ADS_ASSUME_ALIGNED(ptr, CACHELINE_SIZE_BYTES);
        Bitvec bv = createDfuds(container, comp, ptr);
        ptr += bv.template allocatedSizeIn<Group::CacheLine>();
        ADS_ASSUME_ALIGNED(ptr, CACHELINE_SIZE_BYTES);
        return RmmTree(std::move(bv), ptr);
    }

public:
    constexpr static const char name[] = "Succinct Rmq";

    constexpr SuccinctRMQ() = default;

    ADS_CPP20_CONSTEXPR ~SuccinctRMQ() noexcept = default;

    template<typename Container, typename Comp = std::less<typename Container::value_type>>
    ADS_CPP20_CONSTEXPR explicit SuccinctRMQ(const Container& container, Comp comp = Comp())
        : allocation(numAllocatedBytes(container)), rmmTree(createRmmTree(container, comp)), length(container.size()) {
        ADS_ASSUME(allocation.isEnd(rmmTree.rmmArr.end()));
    }

    template<typename T>
    ADS_CPP20_CONSTEXPR SuccinctRMQ(std::unique_ptr<T> ptr, Index length)
        : SuccinctRMQ(Span<const T>(ptr.get(), length)) {}

    template<typename T = Limb>
    ADS_CPP20_CONSTEXPR SuccinctRMQ(std::initializer_list<T> list)
        : SuccinctRMQ(Span<const T>(list.begin(), list.end())) {}

    ADS_CPP20_CONSTEXPR SuccinctRMQ(SuccinctRMQ&&) noexcept = default;

    ADS_CPP20_CONSTEXPR SuccinctRMQ& operator=(SuccinctRMQ&&) noexcept = default;

    ADS_CPP20_CONSTEXPR SuccinctRMQ& operator=(const SuccinctRMQ&) noexcept = delete;

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rmq(Index lower, Index upper) const noexcept {
        assert(lower < upper);
        if (lower + 1 == upper) {
            return lower;
        }
        const Bitvec& dfuds = rmmTree.getBitvector();
        Index x = dfuds.template select<closeParen>(lower);
        Index y = dfuds.template select<closeParen>(upper - 1) + 1;
        ADS_ASSUME(dfuds.getBit(x) == closeParen && dfuds.getBit(y - 1) == closeParen);
        Index w = rmmTree.bitvecRmq(x, y);
        ADS_ASSUME(w >= x && w < y);
        ADS_ASSUME(dfuds.getBit(w) == closeParen);
        // The findOpen call from the paper (https://algo2.iti.kit.edu/download/rmq.pdf, corollary 5.6) is unnecessary:
        // It only handles the case of lca == lower, but in this case, w is equal to x because no descendant of lower can have a smaller excess
        // and lower is the leftmost value with such an excess (see https://doi.org/10.1016/j.jcss.2011.09.002).
        // Therefore, rank<closeParen>(w) already computes the correct answer.
        return dfuds.template rank<closeParen>(w);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index operator()(Index lower, Index upper) const noexcept {
        return rmq(lower, upper);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return length; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index allocatedSizeInBits() const noexcept {
        return Index(rmmTree.size * sizeof(typename RmmTree ::value_type) * 8 + rmmTree.getBitvector().allocatedSizeInBits());
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const RmmTree& getTree() const noexcept { return rmmTree; }
};

} // namespace ads

#endif // ADS_SUCCINCT_RMQ_HPP
