#ifndef ADS_BITVECTOR_HPP
#define ADS_BITVECTOR_HPP

#include "bit.hpp"
#include "bitvec_layouts.hpp"
#include "common.hpp"
#include "rand_access_iter.hpp"

#include <cassert>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <stdexcept>

namespace ads {


template<ADS_LAYOUT_CONCEPT Layout = SimpleLayout<>>
class Bitvector : private Layout {
    using Base = Layout;
    using Base::getBlockCount;
    using Base::getElemRef;
    using Base::getSuperBlockCount;
    using Base::setBlockCount;
    using Base::setSuperBlockCount;

    using typename Base::BlockCount;
    using typename Base::SuperBlockCount;

    const Layout& layout() const noexcept { return static_cast<const Layout&>(*this); }

    Index numBits;

    static Index numElems(Index numBits) noexcept {
        assert(numBits >= 0);
        return roundUpDiv(numBits, 64);
    }

    constexpr static auto getSuperBlockCountFunc
            = [](const auto& bv, Index i) -> SuperBlockCount { return bv.layout().getSuperBlockCount(i); };

    using SuperBlockCountIter = RandAccessIter<Bitvector, decltype(getSuperBlockCountFunc)>;
    SuperBlockCountIter sbIter(Index i) const { return SuperBlockCountIter(*this, getSuperBlockCountFunc, i); }

    Subrange<SuperBlockCountIter> superBlockCountView() const noexcept {
        return Subrange<SuperBlockCountIter>{sbIter(0), sbIter(numSuperBlocks())};
    }


public:
    using Base::allocatedSizeInElems;
    using Base::blockSize;
    using Base::getBit;
    using Base::getElem;
    using Base::numBlocks;
    using Base::numBlocksInSuperBlock;
    using Base::numElems;
    using Base::numElemsInBlock;
    using Base::numElemsInSuperBlock;
    using Base::numSuperBlocks;
    using Base::setBit;
    using Base::setElem;
    using Base::superBlockSize;

    Bitvector() noexcept : numBits(0) {}

    explicit Bitvector(Index numBits) noexcept : Base(numElems(numBits)), numBits(numBits) {}

    Bitvector(Index numBits, Elem fill) : Bitvector(numBits) {
        for (Index i = 0; i < sizeInElems(); ++i) {
            setElem(i, fill);
        }
        buildMetadata();
    }

    explicit Bitvector(Span<Elem> elems) noexcept : Bitvector(elems, elems.size() * 64) {}

    Bitvector(Span<Elem> elems, Index numBits) noexcept : Bitvector(numBits) {
        for (Index i = 0; i < sizeInElems(); ++i) {
            setElem(i, elems[i]);
        }
        buildMetadata();
    }

    explicit Bitvector(std::string_view str, Index base = 2) : Bitvector(Index(str.size()) * log2(Elem(base))) {
        if (base != 2 && base != 4 && base != 16) [[unlikely]] {
            throw std::invalid_argument("base must be one of 2, 4, or 16");
        }
        const Index log2ofBase = log2(Elem(base));
        const std::size_t charsPerElem = 64 / log2ofBase;
        Index i = 0;
        for (Index superblock = 0; superblock < numSuperBlocks(); ++superblock) {
            for (Index j = 0; j < numElemsInSuperBlock(); ++j) {
                if (str.empty()) break;
                std::size_t numToParse = std::min(charsPerElem, str.size());
                std::string_view toParse = str.substr(0, numToParse);
                std::uint64_t res;
                auto err = std::from_chars(toParse.data(), toParse.data() + toParse.size(), res, int(base));
                if (err.ec != std::errc()) [[unlikely]] {
                    throw std::invalid_argument(std::make_error_code(err.ec).message());
                } else if (err.ptr != toParse.data() + toParse.size()) [[unlikely]] {
                    throw std::invalid_argument("invalid character found");
                }
                if (numToParse < charsPerElem) {
                    res <<= (charsPerElem - numToParse) * log2ofBase;
                }
                str.remove_prefix(numToParse);
                res = reverseBits(res);
                setElem(i, res);
                ++i;
            }
            buildRankMetadata(superblock); // TODO: Compute metadata while iterating over elements for cache efficiency
        }
        assert(str.empty());
    }


    void buildRankMetadata(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperBlocks());
        if (superblockIdx == 0) {
            setSuperBlockCount(0, 0);
        }
        auto s = superblockElems(superblockIdx);
        Index inSuperblockSoFar = 0;
        for (Index i = 0; i < s.size(); ++i) {
            if (i % numElemsInBlock() == 0) {
                //                std::cout << "setting block count for block " << i << " to " << inSuperblockSoFar << ", s[i] is "
                //                          << std::hex << s[i] << std::dec << ", popcount " << popcount(s[i]) << std::endl;
                setBlockCount(superblockIdx, i / numElemsInBlock(), inSuperblockSoFar);
            }
            inSuperblockSoFar += popcount(s[i]);
        }
        assert(inSuperblockSoFar <= superBlockSize());
        if (superblockIdx + 1 < numSuperBlocks()) {
            assert(s.size() % numElemsInBlock() == 0);
        }
        // there are numSuperBlocks() + 1 super block counts
        setSuperBlockCount(superblockIdx + 1, getSuperBlockCount(superblockIdx) + inSuperblockSoFar);
    }

    void buildRankMetadata() noexcept {
        for (Index i = 0; i < numSuperBlocks(); ++i) {
            buildRankMetadata(i);
        }
    }

    void buildMetadata() noexcept { buildRankMetadata(); }

    void buildSelectMetadata(Index) noexcept {
        // nothing to do; TODO: Move into Layout class, use CRTP
    }

    [[nodiscard]] Index rankOne(Index pos) const {
        if (pos >= sizeInBits() || pos < 0) [[unlikely]] {
            throw std::invalid_argument("invalid position for rank query");
        }
        return rankOneUnchecked(pos);
    }

    [[nodiscard]] Index numOnes() const noexcept { return getSuperBlockCount(numSuperBlocks()); }

    [[nodiscard]] Index numZeros() const noexcept { return sizeInBits() - numOnes(); }


    /// Unlike rankOne(), this dDoesn't check that `pos` is valid, although that gives close to no measurable performance benefits.
    /// However, the combined ASSUME macros do improve performance by quite a bit (if the compiler couldn't assume that pos >= 0,
    //// performance would actually be significantly lower than with the throwing checks)
    [[nodiscard]] Index rankOneUnchecked(Index pos) const noexcept {
        ADS_ASSUME(0 <= pos);
        ADS_ASSUME(pos < sizeInBits());
        if (pos <= 0) [[unlikely]] {
            return 0;
        }
        --pos;
        ADS_ASSUME(pos >= 0);
        Index elemIdx = pos / 64;
        Elem mask = Elem(-1) >> (63 - pos % 64);
        Index superblockIdx = elemIdx / numElemsInSuperBlock();
        ADS_ASSUME(superblockIdx >= 0);
        Index blockIdx = elemIdx / numElemsInBlock();
        ADS_ASSUME(blockIdx >= 0);
        Index res = getSuperBlockCount(superblockIdx) + getBlockCount(blockIdx);
        ADS_ASSUME(res >= 0);
        for (Index i = blockIdx * numElemsInBlock(); i < elemIdx; ++i) {
            ADS_ASSUME(elemIdx - i < numBlocksInSuperBlock());
            res += popcount(getElem(i));
        }
        return res + popcount(getElem(elemIdx) & mask);
    }

    [[nodiscard]] Index rankZero(Index pos) const { return pos - rankOne(pos); }

    [[nodiscard]] Index rankZeroUnchecked(Index pos) const noexcept { return pos - rankOneUnchecked(pos); }

    template<bool IsOne>
    [[nodiscard]] Index rank(Index pos) const noexcept {
        if constexpr (IsOne) {
            return rankOne(pos);
        } else {
            return rankZero(pos);
        }
    }

private:
    template<bool IsOne>
    [[nodiscard]] Index selectSuperBlockIdx(Index& bitRank) const {
        constexpr Index linearFallbackSize = 8;
        auto rankFunc = [this](Index i) noexcept {
            ADS_ASSUME(i >= 0);
            ADS_ASSUME(i <= numSuperBlocks());
            Index rankOne = getSuperBlockCount(i);
            ADS_ASSUME(rankOne >= 0);
            if constexpr (IsOne) {
                return rankOne;
            } else {
                // if i is numSuperBlocks(), this can be greater than the number of zeros in the bitvector,
                // but that's not a problem
                Index numBitsBefore = i * superBlockSize();
                ADS_ASSUME(rankOne <= numBitsBefore);
                return numBitsBefore - rankOne;
            }
        };
        ADS_ASSUME(bitRank >= 0);
        Index l = bitRank / superBlockSize() + 1; // + 1 because we're searching for the first superblock where the rank is greater than bitRank
        Index u = numSuperBlocks();
        ADS_ASSUME(l > 0);
        ADS_ASSUME(l <= u);
        if (u - l > linearFallbackSize) {
            // set u close to the expected location for iid ones with 50% probability, then increase exponentially
            // until it is an upper bound. Unlike binary search, this starts with a less pessimistic search window and
            // should hopefully be easier on the branch predictor. This improves performance for random values but hurts for especially hard cases.
            u = l;
            do {
                l = u;
                u *= 2;
                if (u >= numSuperBlocks()) {
                    u = numSuperBlocks();
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

    // Idea: For numElemsInBlock >= 4 (and superBlockSize 2^16 bits), the block index can be represented with 1 byte,
    // so using e.g. 32 block indices for selectOne and additional 32 indices for selectZero takes up only 128 Byte = 16
    // Elems, and allows narrowing the search down to 8 blocks on average
    // Also, it's probably safe to optimize for the case of n < 2^48 bits, which means superblock indices can be
    // represented with 32 bit values. Using numSuperBlocks() indices per bit value to select superblocks may be a good
    // idea. All this can be implemented in a layout that derives from the simple layout, in which case selectBlockIdx()
    // etc should be CRTP "virtual" methods of the layout.
    // Even better (?) implementation: For each 2^16 bit superblock, store a 256 bit vector, called a selectBlock, where
    // selectOne(i) gives the block idx of the (i * 256)th one, ie each 256 bit block is represented with one bit.
    // Additionally, for each block where a one is set in the selectBlock, use 8 bit to store the offset between the
    // rank at block start to (i * 256), which allows figuring out which bit in the block was the (i * 256)th bit in the
    // superblock. Do the same for selectZero but reuse the same 256 Byte array to store block offsets for select1 and
    // select0, which works because the sum of set bits in the select1 and select0 selectBlocks is at most 256.
    // Then, it's not even necessary to store block counts, although that would make select be linear within a
    // superblock in the worst case (for uniformly iid ones with probability p, the expected number of blocks to search
    // is min(1 / p, 2^8). It's also possible to explicitly store the index of all set bits using 16 bit per entry in
    // the same array if the number of set bits is at most 128 (or, equivalently, the number of unset bits).
    // However, this doesn't help if the superblock consists of 128 ones, then only zeros until the last bit, which is
    // a 1. Even better idea: If the selectBlock bit is a 0, store the rank % 256 instead, which allows binary search in
    // such cases. This effectively doubles the required memory, totalling to 1/32th + 1/128th of the size in bits.
    // Also, this can be applied to the entire bitvector, not just individual superblocks. The recursive bit vector can
    // be stored in a more space consuming manner.
    // TODO: Also, it may be worth investigating if binary search works better if the range isn't split in two equal
    // sized halves but instead the requested bitRank as a fraction of the total rank in the (bv|superblock|block) is
    // used to generate the pivot element. The problem of this approach is that it can be very inefficient for
    // adversarial patters such as a superblock with 100 ones, all of which are in the last 100 bits and a selectOne(0)
    // query. Maybe some combination with splitting in equal-sized halves can give reasonable worst-case guarantees?
    template<bool IsOne>
    [[nodiscard]] Index selectBlockIdx(Index& bitRank, Index superBlockIdx) const noexcept {
        auto rankFunc = [this](Index i) noexcept {
            Index rankOne = getBlockCount(i);
            if constexpr (IsOne) {
                return rankOne;
            } else {
                Index numBitsBefore = (i % numBlocksInSuperBlock()) * blockSize();
                ADS_ASSUME(rankOne <= numBitsBefore);
                return numBitsBefore - rankOne;
            }
        };
        constexpr Index linearFallbackSize = 2 * sizeof(Elem) / sizeof(BlockCount);
        ADS_ASSUME(superBlockIdx >= 0);
        ADS_ASSUME(superBlockIdx < numSuperBlocks());
        ADS_ASSUME(bitRank >= 0);
        // we're searching for the first block with count strictly greater than bitRank
        Index l = superBlockIdx * numBlocksInSuperBlock() + 1;
        Index u = std::min(l + numBlocksInSuperBlock() - 1, numBlocks());
        ADS_ASSUME(u >= l);
        ADS_ASSUME(u - l < numBlocksInSuperBlock());
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
        ADS_ASSUME(l > superBlockIdx * numBlocksInSuperBlock());
        for (Index i = l; i < u; ++i) {
            ADS_ASSUME(i == 0 || getBlockCount(i) >= getBlockCount(i - 1));
            if (rankFunc(i) > bitRank) {
                bitRank -= rankFunc(i - 1);
                return i - 1;
            }
        }
        bitRank -= rankFunc(u - 1);
        return u - 1;
    }

    template<bool IsOne>
    [[nodiscard]] Index selectElemIdx(Index& bitRank, Index blockIdx) const noexcept {
        Index first = blockIdx * numElemsInBlock();
        auto rankFunc = [this](Index i) noexcept {
            Elem e = getElem(i);
            if constexpr (IsOne) {
                return popcount(e);
            } else {
                return 64 - popcount(e);
            }
        };
        for (Index i = first; i < numElems(); ++i) {
            Index rank = rankFunc(i);
            ADS_ASSUME(i - first < numElemsInBlock());
            ADS_ASSUME(rank >= 0);
            if (rank > bitRank) {
                return i;
            }
            bitRank -= rank;
            ADS_ASSUME(bitRank >= 0);
        }
        return numElems() - 1;
    }

    template<bool IsOne>
    [[nodiscard]] Index selectBitIdx(Elem elem, Index bitIndex) const noexcept {
        if constexpr (IsOne) {
            return elemSelect(elem, bitIndex);
        } else {
            return elemSelect(~elem, bitIndex);
        }
    }

public:
    template<bool IsOne>
    [[nodiscard]] Index select(Index bitRank) const {
        if (bitRank < 0 || bitRank >= sizeInBits()) [[unlikely]] {
            throw std::invalid_argument("invalid rank for select query: " + std::to_string(bitRank));
        }
        //        Index removeThis = bitRank;
        Index superBlockIdx = selectSuperBlockIdx<IsOne>(bitRank);
        Index blockIdx = selectBlockIdx<IsOne>(bitRank, superBlockIdx);
        Index elemIdx = selectElemIdx<IsOne>(bitRank, blockIdx);
        Index bitIdx = selectBitIdx<IsOne>(getElem(elemIdx), bitRank);
        return elemIdx * 64 + bitIdx;
        //        bitRank = removeThis;
        //        Index lower = elemIdx * 64;
        //        Index upper = std::min(lower + 64, sizeInBits());
        //        ADS_ASSUME(lower < upper);
        //        //        Index lower = 0, upper = sizeInBits(), mid = -1;
        //        while (upper - lower > 1) { // TODO: linear fallback
        //            Index mid = (lower + upper) / 2;
        //            Index r = rank<IsOne>(mid);
        //            if (r <= bitRank) {
        //                lower = mid;
        //            } else {
        //                upper = mid;
        //            }
        //        }
        //        return rank<IsOne>(lower) == bitRank && getBit(lower) == IsOne ? lower : -1;
    }

    /// Return the position `i` such that `getBit(i)` returns the `rank`th one, counting from 0.
    [[nodiscard]] Index selectOne(Index rank) const { return select<true>(rank); }

    /// Return the position `i` such that `getBit(i)` returns the `rank`th zero, counting from 0
    [[nodiscard]] Index selectZero(Index rank) const { return select<false>(rank); }


    // Idea: superblock size is 4 Elems = 32 Byte = 256 bit, use 1 Byte to store number of 0s in superblock -- numbers
    // in range [0, 8 * 24] with s' = 256, we have s = 2^4 = 16 = (log n) / 2, therefore n = 2 ^ 32 n + 8 + n / 8 <= 64
    // <=> n <= 56 * 8 / 9 => 48, ergo 48 + 8 + 6 bytes, 48 Byte = cache efficient layout: superblock bitvector (32
    // Byte) + Superblock count (8 Byte) + block count (1 Byte) * 4 == 44 Byte, wasting 20 Byte or : Superblock size is
    // less than 2^16 bit = 65636 bit = 4096 Byte = 512 Elems because s = (log n) / 2 <= 32, s' = s * s <= 1024

    [[nodiscard]] Index sizeInBits() const noexcept { return numBits; }

    [[nodiscard]] Index sizeInElems() const noexcept {
        assert(numElems() == roundUpDiv(numBits, 64));
        return numElems();
    }

    [[nodiscard]] Span<Elem> superblockElems(Index superblockIdx) noexcept {
        Index size = superblockIdx == numSuperBlocks() - 1 ? (numElems() - 1) % numElemsInSuperBlock() + 1 :
                                                             numElemsInSuperBlock();
        return Span<Elem>(&getElemRef(superblockIdx * numElemsInSuperBlock()), size);
    }
    [[nodiscard]] Span<const Elem> superblockElems(Index superblockIdx) const noexcept {
        Index size = superblockIdx == numSuperBlocks() - 1 ? (numElems() - 1) % numElemsInSuperBlock() + 1 :
                                                             numElemsInSuperBlock();
        return Span<const Elem>(&getElemRef(superblockIdx * numElemsInSuperBlock()), size);
    }

    [[nodiscard]] Index numAllocatedBits() const noexcept { return allocatedSizeInElems() * sizeof(Elem) * 8; }

    constexpr static auto getBitFunc = [](const auto& bv, Index i) -> bool { return bv.getBit(i); };
    constexpr static auto getElemFunc = [](const auto& bv, Index i) -> Elem { return bv.getElem(i); };

    using BitIter = RandAccessIter<Bitvector, decltype(getBitFunc)>;

    BitIter bitIter(Index i) const { return BitIter(*this, getBitFunc, i); }

    Subrange<BitIter> bitView() const noexcept { return Subrange<BitIter>{bitIter(0), bitIter(numBits)}; }

    using ElemIter = RandAccessIter<Bitvector, decltype(getElemFunc)>;

    ElemIter elemIter(Index i) const { return ElemIter(*this, getElemFunc, i); }

    Subrange<ElemIter> elemView() const noexcept { return Subrange<ElemIter>{elemIter(0), elemIter(numElems())}; }

    friend std::ostream& operator<<(std::ostream& os, const Bitvector& bv) {
        std::copy(bv.bitView().begin(), bv.bitView().end(), std::ostream_iterator<bool>(os, ""));
        return os;
    }

    [[nodiscard]] std::string toString() const noexcept {
        std::string res;
        for (bool b : bitView()) {
            if (b) {
                res += '1';
            } else {
                res += '0';
            }
        }
        return res;
    }
};


template<ADS_LAYOUT_CONCEPT L1, ADS_LAYOUT_CONCEPT L2>
bool operator==(const Bitvector<L1>& lhs, const Bitvector<L2>& rhs) {
    return lhs.sizeInBits() == rhs.sizeInBits()
           && std::equal(lhs.elemView().begin(), lhs.elemView().end(), rhs.elemView().begin());
}
template<ADS_LAYOUT_CONCEPT L1, ADS_LAYOUT_CONCEPT L2>
bool operator!=(const Bitvector<L1>& lhs, const Bitvector<L2>& rhs) {
    return !(rhs == lhs);
}

#ifdef ADS_HAS_CPP20

template<ADS_LAYOUT_CONCEPT L1, ADS_LAYOUT_CONCEPT L2>
std::strong_ordering operator<=>(const Bitvector<L1>& lhs, const Bitvector<L2>& rhs) noexcept {
    if (lhs.sizeInBits() != rhs.sizeInBits()) {
        return std::lexicographical_compare_three_way(
                lhs.bitView().begin(), lhs.bitView().end(), rhs.bitView().begin(), rhs.bitView().end());
    }
    return std::lexicographical_compare_three_way(
            lhs.elemView().begin(), lhs.elemView().end(), rhs.elemView().begin(), rhs.elemView().end());
}

#else

template<ADS_LAYOUT_CONCEPT L1, ADS_LAYOUT_CONCEPT L2>
bool operator<(const Bitvector<L1>& lhs, const Bitvector<L2>& rhs) {
    if (lhs.sizeInBits() != rhs.sizeInBits()) {
        return std::lexicographical_compare(
                lhs.bitView().begin(), lhs.bitView().end(), rhs.bitView().begin(), rhs.bitView().end());
    }
    return std::lexicographical_compare(
            lhs.elemView().begin(), lhs.elemView().end(), rhs.elemView().begin(), rhs.elemView().end());
}
template<ADS_LAYOUT_CONCEPT L1, ADS_LAYOUT_CONCEPT L2>
bool operator>(const Bitvector<L1>& lhs, const Bitvector<L2>& rhs) {
    return rhs < lhs;
}
template<ADS_LAYOUT_CONCEPT L1, ADS_LAYOUT_CONCEPT L2>
bool operator<=(const Bitvector<L1>& lhs, const Bitvector<L2>& rhs) {
    return !(rhs < lhs);
}
template<ADS_LAYOUT_CONCEPT L1, ADS_LAYOUT_CONCEPT L2>
bool operator>=(const Bitvector<L1>& lhs, const Bitvector<L2>& rhs) {
    return !(lhs < rhs);
}

#endif

} // namespace ads

#endif // ADS_BITVECTOR_HPP
