
#ifndef BITVECTOR_SUCCINCT_RMQ_HPP
#define BITVECTOR_SUCCINCT_RMQ_HPP

#include "bitvector.hpp"
#include "common.hpp"
#include <algorithm>
#include <stack>
namespace ads {

/// This file implements the succinct rmq data structure from https://algo2.iti.kit.edu/download/rmq.pdf,
/// using the range min-max tree from https://arxiv.org/pdf/1601.06939.pdf for solving the excess rmq problem on bitvectors.

constexpr static bool openParen = true;
constexpr static bool closeParen = false;

template<typename T = Elem, Index BlockSize = 512>// TODO: Make sure the bitvector allocates cacheline-aligned
struct RangeMinMaxTree {
    static_assert(BlockSize <= std::numeric_limits<T>::max() / 2);
    static_assert(BlockSize % 64 == 0);

    static constexpr Index blockSize = BlockSize;

    std::unique_ptr<T[]> rmmArr = nullptr;
    Bitvector<> bv;
    Index numLeaves = 0;
    Index size = 0;

    RangeMinMaxTree() = default;

    [[nodiscard]] Index leafIdxInArr(Index leafNum) const noexcept {
        return size - numLeaves + leafNum;
    }

    /// The reverse of leafIdxInArr
    [[nodiscard]] Index leafNum(Index leafIdxInArr) const noexcept {
        assert(leafIdxInArr >= size - numLeaves && leafIdxInArr < size);
        return leafIdxInArr - (size - numLeaves);
    }

    [[nodiscard]] bool isRightChild(Index node) const noexcept {
        return node % 2 == 1;
    }

    // TODO: Different implementation: Store for each leaf the position of its minimum and the minimum excess in block (>= -2^BlockLength, <= 2^BlockLength),
    // which can be done with 16 bits per leaf if the block size is no greater than 256 (min. excess in block can be relative to next block to keep numbers small).
    // Implement and measure if faster (for the same block size, this additionally uses 16 bits per leaf, which is roughly 8 bits per node but may cause cache misses)

    explicit RangeMinMaxTree(Bitvector<>&& bitvector) noexcept : bv(std::move(bitvector)) {
        numLeaves = roundUpDiv(bv.sizeInBits(), BlockSize);
        assert(numLeaves > 0);
        size = numLeaves + (Elem(1) << roundUpLog2(Elem(numLeaves)));// the missing -1 makes this a 1-indexed heap
        assert(leafIdxInArr(0) % 2 == 0 || leafIdxInArr(0) == 1);
        rmmArr = makeUniqueForOverwrite<T>(size);
        Index minExcessInBlock = 1;
        Index excess = 0;
        for (Index i = 0; i < bv.sizeInBits(); ++i) {
            if (i > 0 && i % BlockSize == 0) {
                Index l = leafIdxInArr(i / BlockSize - 1);
                rmmArr[l] = minExcessInBlock;
                while (isRightChild(l)) {
                    rmmArr[l / 2] = std::min(rmmArr[l - 1], rmmArr[l]);
                    l /= 2;
                }
                minExcessInBlock = excess + 1;
            }
            bool bit = bv.getBit(i);
            if (bit == openParen) {// TODO: There's probably a better way to do this than to read each bit individually - read bytes and use lookup table?
                ++excess;
            } else {
                --excess;
                minExcessInBlock = std::min(minExcessInBlock, excess);
            }
        }
        assert(minExcessInBlock == 0);
        Index v = leafIdxInArr(numLeaves - 1);
        rmmArr[v] = 0;
        for (; v > 0; v /= 2) {
            if (isRightChild(v)) {
                rmmArr[v / 2] = std::min(rmmArr[v], rmmArr[v - 1]);
            } else {
                rmmArr[v / 2] = rmmArr[v];
            }
        }
    }
    void printTree() const {
        std::cout << " tree stored in array with " << size << " values, " << numLeaves << " leaves:" << std::endl;
        for (Index i = 1; i < this->size; ++i) {
            if (i == this->leafIdxInArr(0)) { std::cout << "(leafs) "; }
            std::cout << this->rmmArr[i] << " ";
            if (this->isRightChild(i)) { std::cout << ' '; }
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

    MinRes findMinInBlockImpl(Index lower, Index upper, Index& excessSoFar) const noexcept {
        if (lower == upper) { return noRes; }
        Index excess = excessSoFar;
        Index minExcess = excess + 1;
        Index minPos = lower;
        for (Index i = lower; i < upper; ++i) {
            if (bv.getBit(i) == openParen) {
                ++excess;
            } else if (--excess < minExcess) {// choose leftmost index
                minExcess = excess;
                minPos = i;
            }
        }
        excessSoFar = excess;
        return {minPos, minExcess};
    }

    MinRes findMinInBlock(Index lower, Index upper, Index& excessUntilUpper, Index& globalOffset) const noexcept {
        Index blockBegin = lower / BlockSize * BlockSize;
        Index blockEnd = blockBegin + BlockSize;
        blockEnd = std::min(blockEnd, bv.sizeInBits());
        upper = std::min(upper, bv.sizeInBits());
        assert(lower <= upper);
        if (upper == lower) {
            Index ignored = 0;
            Index minExcessInBlock = findMinInBlockImpl(blockBegin, blockEnd, ignored).minExcess;
            globalOffset = rmmArr[leafIdxInArr(lower / BlockSize)] - minExcessInBlock;
            assert(globalOffset >= 0 || blockBegin == blockEnd);
            return noRes;
        }
        Index inBlockExcess = 0;
        Index minExcessInBlock = findMinInBlockImpl(blockBegin, lower, inBlockExcess).minExcess;
        MinRes res = findMinInBlockImpl(lower, upper, inBlockExcess);
        excessUntilUpper = inBlockExcess;
        assert(res.pos >= lower && res.pos < upper);
        minExcessInBlock = std::min(minExcessInBlock, res.minExcess);
        Index minExcessAfterQuery = findMinInBlockImpl(upper, blockEnd, inBlockExcess).minExcess;
        minExcessInBlock = std::min(minExcessInBlock, minExcessAfterQuery);
        globalOffset = rmmArr[leafIdxInArr(lower / BlockSize)] - minExcessInBlock;
        return {res.pos, res.minExcess + globalOffset};
    }

    // Idea: find min extent within block and use that together with the heap's stored value of the global min extent in this block
    // to figure out the extent at the start of the block, which allows calculating the global min extent within the query range
    [[nodiscard]] MinRes findMinInBlock(Index lower, Index upper) const noexcept {
        Index ignored = 0;
        return findMinInBlock(lower, upper, ignored, ignored);
    }


    /// \brief Find the minimum element in the half-open range of blocks [\p lowerBlockIdx, \p upperBlockIdx)
    /// \param lowerBlockIdx
    /// \param upperBlockIdx
    /// \return
    [[nodiscard]] MinRes findMinInTree(Index lowerBlockIdx, Index upperBlockIdx) const noexcept {
        assert(lowerBlockIdx <= upperBlockIdx + 1);
        if (upperBlockIdx <= lowerBlockIdx) {
            return noRes;
        } else if (lowerBlockIdx + 1 == upperBlockIdx) {
            return findMinInBlock(lowerBlockIdx * BlockSize, upperBlockIdx * BlockSize);
        }
        Index lower = leafIdxInArr(lowerBlockIdx);
        Index upper = leafIdxInArr(upperBlockIdx - 1);
        assert(log2(Elem(lower)) == log2(Elem(upper)));
        Index h = log2(Elem(lower) ^ Elem(upper)) + 1;
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
        Index minPos = rmmArr[rightMinPos] < rmmArr[leftMinPos] ? rightMinPos : leftMinPos;// choose left on ties
        while (minPos < lower) {
            minPos *= 2;
            assert(minPos <= upper);
            if (rmmArr[minPos + 1] < rmmArr[minPos]) {
                ++minPos;
            }
            assert(minPos <= upper);
        }
        Index i = leafNum(minPos);
        assert(std::min_element(rmmArr.get() + lower, rmmArr.get() + upper + 1) - rmmArr.get() == minPos);
        return findMinInBlock(i * BlockSize, (i + 1) * BlockSize);
    }

    [[nodiscard]] MinRes rmqImpl(Index lower, Index upper) const noexcept {
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
        MinRes lowerOrTreeCandidate = treeCandidate.minExcess < lowerCandidate.minExcess ? treeCandidate : lowerCandidate;// take lower when tied
        MinRes upperCandidate = findMinInBlock(prevUpperBlockIdx * BlockSize, upper);
        assert((upperCandidate.pos >= lower && upperCandidate.pos < upper) || upperCandidate.pos == noRes.pos);
        return upperCandidate.minExcess < lowerOrTreeCandidate.minExcess ? upperCandidate : lowerOrTreeCandidate;// don't take upper when tied
    }

    [[nodiscard]] Index bitvecRmq(Index lower, Index upper) const noexcept {
        return rmqImpl(lower, upper).pos;
    }

    [[nodiscard]] Index findLastWithExcess(Index last, Index excess) const {
        assert(excess >= 0);
        Index i = last;
        if (excess == 0) {
            assert(bv.getBit(last) == openParen);
            return last;
        }
        while (true) {
            --i;
            assert(i >= (last - 1) / BlockSize * BlockSize);
            if (bv.getBit(i) == closeParen) {
                ++excess;
            } else if (--excess <= 0) {
                break;
            }
        }
        return i;
    }

    // The rightmost index where the excess is <= excess must be achieved by a closing parenthesis (unless it's the first bit in the bitvector).
    // Find that parenthesis and add 1 to find the matching open parenthesis.
    [[nodiscard]] Index findOpen(Index i) const noexcept {
        if (i == bv.sizeInBits() - 1) { return 0; }// don't try to find a bit with excess -1
        Index excess = 0;
        Index excessOfBlockBegin = 0;
        assert(bv.getBit(i) == closeParen);
        Index blockBegin = i / BlockSize * BlockSize;
        MinRes minInBlockBeforeI = findMinInBlock(blockBegin, i, excess, excessOfBlockBegin);
        excess += excessOfBlockBegin - 1;// -1 because the closing ')' must also be counted
        assert(excess >= Index(rmmArr[leafIdxInArr(i / BlockSize)]));
        assert(excess >= 0);
        assert(minInBlockBeforeI.minExcess <= excess + 1 || minInBlockBeforeI.minExcess == noRes.minExcess);
        if (minInBlockBeforeI.minExcess <= excess) {
            assert(minInBlockBeforeI.pos < i && minInBlockBeforeI.pos >= blockBegin);
            return findLastWithExcess(i, 1);
        }
        // move up in tree
        assert(i >= BlockSize);
        Index v = leafIdxInArr(i / BlockSize);
        while (true) {
            if (isRightChild(v) && Index(rmmArr[v - 1]) <= excess) {
                --v;
                break;
            }
            v /= 2;
            assert(v > 1);
        }
        assert(v >= 1);
        assert(v == 1 || !isRightChild(v));
        // move down in tree
        Index firstLeaf = leafIdxInArr(0);
        while (v < firstLeaf) {
            assert(Index(rmmArr[v]) <= excess);
            v *= 2;
            if (Index(rmmArr[v + 1]) <= excess) {
                ++v;
            }
        }
        assert(Index(rmmArr[v]) <= excess);
        // find the last position with the required excess in the current block
        assert(leafNum(v) < i / BlockSize);
        Index blockEnd = (leafNum(v) + 1) * BlockSize;
        Index excessInBlock = 0;
        findMinInBlock(blockEnd - BlockSize, blockEnd, excessInBlock, excessOfBlockBegin);
        Index currentExcess = excessInBlock + excessOfBlockBegin;
        assert(currentExcess >= Index(rmmArr[v]));
        assert(blockBegin == blockEnd || Index(rmmArr[v + 1]) > excess);
        assert(currentExcess >= excess);
        return findLastWithExcess(blockEnd, currentExcess - excess);
    }

    [[nodiscard]] const Bitvector<>& getBitvector() const noexcept { return bv; }
};


class SuccinctRMQ {

    RangeMinMaxTree<> rmmTree = RangeMinMaxTree<>();
    Index length = 0;

public:
    constexpr static const char name[] = "Succinct Rmq";

    SuccinctRMQ() = default;

    template<typename Container, typename Comp = std::less<typename Container::value_type>>
    explicit SuccinctRMQ(const Container& container, Comp comp = Comp()) : length(container.size()) {
        Bitvector<> dfuds(container.size() * 2 + 2);
        Index bvIdx = dfuds.sizeInBits();
        std::stack<Index, std::vector<Index>> stack;
        for (Index i = container.size() - 1; i + 1 > 0; --i) {
            assert(bvIdx > 1);
            dfuds.setBit(--bvIdx, closeParen);
            while (!stack.empty() && comp(container[i], container[stack.top()])) {
                assert(bvIdx > 1);
                dfuds.setBit(--bvIdx, openParen);
                stack.pop();
            }
            stack.push(i);
        }
        dfuds.setBit(--bvIdx, closeParen);
        while (!stack.empty()) {
            stack.pop();
            dfuds.setBit(--bvIdx, openParen);
        }
        assert(bvIdx == 1);
        dfuds.setBit(0, openParen);
        for (Index i = 0; i < dfuds.numSuperblocks(); ++i) {
            dfuds.buildRankMetadata(i);
        }
        rmmTree = RangeMinMaxTree<>(std::move(dfuds));
    }
    template<typename T>
    SuccinctRMQ(std::unique_ptr<T> ptr, Index length) : SuccinctRMQ(Span<const T>(ptr.get(), length)) {}
    template<typename T = Elem>
    SuccinctRMQ(std::initializer_list<T> list) : SuccinctRMQ(Span<const T>(list.begin(), list.end())) {}

    [[nodiscard]] Index rmq(Index lower, Index upper) const noexcept {
        assert(lower < upper);
        if (lower + 1 == upper) { return lower; }
        const Bitvector<>& dfuds = rmmTree.getBitvector();
        Index x = dfuds.select<closeParen>(lower);
        Index y = dfuds.select<closeParen>(upper - 1) + 1;
        assert(dfuds.getBit(x) == closeParen && dfuds.getBit(y - 1) == closeParen);
        Index w = rmmTree.bitvecRmq(x, y);
        assert(w >= x && w < y);
        assert(dfuds.getBit(w) == closeParen);
        if (dfuds.rank<closeParen>(rmmTree.findOpen(w)) == lower + 1) {
            return lower;
        }
        return dfuds.rank<closeParen>(w);
    }

    [[nodiscard]] Index operator()(Index lower, Index upper) const noexcept {
        return rmq(lower, upper);
    }

    [[nodiscard]] Index size() const noexcept {
        return length;
    }

    [[nodiscard]] const RangeMinMaxTree<>& getTree() const noexcept {
        return rmmTree;
    }
};

}// namespace ads

#endif//BITVECTOR_SUCCINCT_RMQ_HPP
