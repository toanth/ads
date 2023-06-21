#ifndef ADS_BITVECTOR_HPP
#define ADS_BITVECTOR_HPP

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
        if (superblockIdx < numSuperBlocks() - 1) {
            assert(s.size() % numElemsInBlock() == 0);
            setSuperBlockCount(superblockIdx + 1, getSuperBlockCount(superblockIdx) + inSuperblockSoFar);
        }
    }

    void buildRankMetadata() noexcept {
        for (Index i = 0; i < numSuperBlocks(); ++i) {
            buildRankMetadata(i);
        }
    }

    void buildMetadata() noexcept { buildRankMetadata(); }

    void buildSelectMetadata(Index) noexcept {
        // idea: store array of superblock start indices
        // superblock size 1024 zeros, minimum size 1024 bit = 2^10 bit = 128 Byte = 16 Elem,
        // maximum size 65536 bit = 2^16 bit = 8192 Byte = 1024 Elems; memory usage <= 8 Byte = 64 bit per 2^10 bits =
        // 1/16th of bv block size 2^7 = 128 zeros; 2^3 = 8 blocks in a superblock; memory usage 32 bit per block => 256
        // bits = 32 Byte = 4 Elem per superblock maximum block size = roughly 2^16 bits
        // -- idea: only store elem idx, which saves log2(64) = 6 bits, use those to store number of zeros in Elem
        // before the original position for blocks <= 128 * 16 * 8 bits = 16384 bits = 2048 Bytes = 256 Elem, don't
        // store anything and compute answer by looking at the bv, possibly use rank queries else, store all answers
        // naively
        // -- what about select 0? reuse select 1? how?
    }

    [[nodiscard]] Index rankOne(Index pos) const {
        if (pos >= sizeInBits() || pos < 0) [[unlikely]] {
            throw std::invalid_argument("invalid position for rank query");
        }
        return rankOneUnchecked(pos);
    }

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

    template<bool IsOne>
    [[nodiscard]] Index selectSuperBlockIdx(Index bitRank) const {
        constexpr Index linearFallbackSize = 64;
        auto rankFunc = [this](Index i) noexcept {
            ADS_ASSUME(i >= 0);
            ADS_ASSUME(i < numSuperBlocks());
            Index rankOne = getSuperBlockCount(i);
            ADS_ASSUME(rankOne >= 0);
            if constexpr (IsOne) {
                return rankOne;
            } else {
                Index numBitsBefore = i * superBlockSize();
                ADS_ASSUME(rankOne <= numBitsBefore);
                return numBitsBefore - rankOne;
            }
        };
        ADS_ASSUME(bitRank >= 0);
        Index l = bitRank / superBlockSize();
        Index u = numSuperBlocks();
        ADS_ASSUME(l <= u);
        if (u - l > linearFallbackSize) {
            // set u close to the expected location for iid ones with 50% probability, then increase exponentially
            // until it is an upper bound. Unlike binary search, this starts with a less pessimistic search window and
            // should hopefully be easier on the branch predictor
            u = l + 1; // add one in case l == 0
            do {
                u *= 2;
                if (u >= numSuperBlocks()) {
                    u = numSuperBlocks();
                    break;
                }
            } while (rankFunc(u) <= bitRank);
            while (l + linearFallbackSize < u) {
                ADS_ASSUME(0 <= l);
                ADS_ASSUME(l < u);
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
        for (Index i = l; i < u; ++i) {
            if (rankFunc(i) > bitRank) {
                return i - 1;
            }
        }
        return u - 1;
    }

    template<bool IsOne>
    [[nodiscard]] Index select(Index bitRank) const {
        if (bitRank < 0 || bitRank >= sizeInBits()) [[unlikely]] {
            throw std::invalid_argument("invalid rank for select query: " + std::to_string(bitRank));
        }
        Index superBlockIdx = selectSuperBlockIdx<IsOne>(bitRank);
        Index lower = superBlockIdx * superBlockSize();
        Index upper = std::min(lower + superBlockSize(), sizeInBits());
        //        Index lower = 0, upper = sizeInBits(), mid = -1;
        while (upper - lower > 1) { // TODO: linear fallback
            Index mid = (lower + upper) / 2;
            Index r = rank<IsOne>(mid);
            if (r <= bitRank) {
                lower = mid;
            } else {
                upper = mid;
            }
        }
        return rank<IsOne>(lower) == bitRank && getBit(lower) == IsOne ? lower : -1;
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
