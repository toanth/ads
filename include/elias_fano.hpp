#ifndef BITVECTOR_ELIAS_FANO_HPP
#define BITVECTOR_ELIAS_FANO_HPP

#include "bitvector.hpp"
#include <cmath>
namespace ads {

template<typename Number = std::uint64_t, ADS_LAYOUT_CONCEPT BitvecLayout = SimpleLayout<>>
class EliasFano {

    Index numInts = 0;
    Index numLowerBitsPerNumber = 0;
    Index lowerBitMask;
    Index completeSizeInBits;
    Bitvector<BitvecLayout> upper;
    BitwiseAccess<dynSize> lower;

private:
    static constexpr Index numBitsInNumber = 8 * sizeof(Number);

    struct CreateWithSizeTag {};

    EliasFano(Index numInts, CreateWithSizeTag)
        : numInts(numInts), numLowerBitsPerNumber(numBitsInNumber - Index(std::ceil(std::log2(numInts)))) {
        Index lowerSizeInElems = roundUpDiv(numLowerBitsPerNumber * numInts, 8 * sizeof(Elem));
        lower.vec = makeUniqueForOverwrite<Elem>(lowerSizeInElems);
        lowerBitMask = (Number(1) << numLowerBitsPerNumber) - 1;
        upper = Bitvector<BitvecLayout>(numInts + (1 << getNumUpperBitsPerNumber()) + 1);
        completeSizeInBits = 8 * (sizeof(*this) + (lowerSizeInElems + upper.completeSizeInElems()) * sizeof(Elem));
        // max: maxIdx = n - 1, max value = 2^numbits - 1
    }


public:
    template<typename Integer, typename = std::enable_if_t<std::is_convertible_v<Integer, Index>>>
    explicit EliasFano(Integer numInts) : EliasFano(numInts, CreateWithSizeTag{}) {
    }

    template<typename Range, typename = std::void_t<decltype(maybe_ranges::begin(std::declval<Range&>()))>>// no concepts in C++17
    explicit EliasFano(const Range& numbers)
        : EliasFano(maybe_ranges::begin(numbers), maybe_ranges::end(numbers)) {
    }

    EliasFano(std::initializer_list<Number> list) : EliasFano(maybe_ranges::begin(list), maybe_ranges::end(list)) {
    }

    // TODO: Instead of the worst case assumpption that max(first, last) == (1 << numBits) - 1, actually check the last element
    template<typename ForwardIter, typename Sentinel>
    EliasFano(ForwardIter first, Sentinel last) noexcept : EliasFano(maybe_ranges::distance(first, last), CreateWithSizeTag{}) {
        Index lowerBitIdx = 0;
        Elem currentUpperEntry = 0;
        Index lastUpperElemIdx = 0;
        Index elemIdx = 1;
        while (first != last) {
            Elem upperPart = Elem(*first) >> numLowerBitsPerNumber;
            Elem newBitIdx = upperPart + elemIdx;
            Index currentUpperElemIdx = newBitIdx / 64;
            if (currentUpperElemIdx > lastUpperElemIdx) {
                upper.setElem(lastUpperElemIdx, currentUpperEntry);
                for (Index k = lastUpperElemIdx + 1; k < currentUpperElemIdx; ++k) {
                    upper.setElem(k, 0);
                }
                lastUpperElemIdx = currentUpperElemIdx;
                currentUpperEntry = 0;
            }
            currentUpperEntry |= Elem(1) << (63 - (newBitIdx % 64));// the order of bits in an element is reversed
            lower.setBits(lowerBitIdx, Elem(*first), numLowerBitsPerNumber);
            lowerBitIdx += numLowerBitsPerNumber;
            ++elemIdx;
            ++first;
        }
        upper.setElem(lastUpperElemIdx, currentUpperEntry);
        for (++lastUpperElemIdx; lastUpperElemIdx < upper.sizeInElems(); ++lastUpperElemIdx) {
            upper.setElem(lastUpperElemIdx, 0);
        }
        for (Index i = 0; i < upper.numSuperblocks(); ++i) {
            upper.buildRankMetadata(i);
        }
        upper.buildSelectMetadata(-1);
    }

    [[nodiscard]] const Bitvector<BitvecLayout>& getUpper() const noexcept {
        return upper;
    }

    [[nodiscard]] const BitwiseAccess<dynSize>& getLower() const noexcept {
        return lower;
    }

    [[nodiscard]] Elem getUpperPart(Index i) const {
        return upper.selectOne(i) - i - 1;
    }

    [[nodiscard]] Elem getLowerPart(Index i) const {
        return lower.getBits(i * numLowerBitsPerNumber, numLowerBitsPerNumber);
    }

    [[nodiscard]] Number get(Index i) const {
        if (i < 0 || i >= size()) {
            throw std::invalid_argument("EliasFano::get() Index out of range");
        }
        Elem upperPart = getUpperPart(i);
        Elem lowerPart = getLowerPart(i);
        return Number(lowerPart + (upperPart << numLowerBitsPerNumber));
    }

    [[nodiscard]] Number operator[](Index i) const {
        return get(i);
    }

    [[nodiscard]] Number predecessor(Number n) const {
        Number upperSearchBits = n >> numLowerBitsPerNumber;
        Number lowerSearchBits = n & lowerBitMask;
        Index first = upper.selectZero(upperSearchBits) - upperSearchBits;
        Index last = upper.selectZero(upperSearchBits + 1) - upperSearchBits - 1;
        for (Index i = last - 1; i >= first; --i) {
            Number lowerBits = lower.getBits(i * numLowerBitsPerNumber, numLowerBitsPerNumber);
            if (lowerBits <= lowerSearchBits) {
                return lowerBits + (upperSearchBits << numLowerBitsPerNumber);
            }
        }
        // there was no element with the same upper part, so return the greatest element with a smaller upper part
        if (first == 0) [[unlikely]] {
            throw std::invalid_argument("no predecessor found for " + std::to_string(n));
        }
        return get(first - 1);
    }

    [[nodiscard]] Number successor(Number n) const {
        Number upperSearchBits = n >> numLowerBitsPerNumber;
        Number lowerSearchBits = n & lowerBitMask;
        Index first = upper.selectZero(upperSearchBits) - upperSearchBits;
        Index last = upper.selectZero(upperSearchBits + 1) - upperSearchBits - 1;
        for (Index i = first; i != last; ++i) {
            Number lowerBits = lower.getBits(i * numLowerBitsPerNumber, numLowerBitsPerNumber);
            if (lowerBits >= lowerSearchBits) {
                return lowerBits + (upperSearchBits << numLowerBitsPerNumber);
            }
        }
        // there was no element with the same upper part, so return the greatest element with a smaller upper part
        if (last == numInts) [[unlikely]] {
            throw std::invalid_argument("no successor found for " + std::to_string(n));
        }
        return get(last);
    }

    [[nodiscard]] Index size() const noexcept {
        return numInts;
    }

    [[nodiscard]] Index spaceInBits() const noexcept {
        return completeSizeInBits;
    }

    [[nodiscard]] Index getNumLowerBitsPerNumber() const noexcept {
        return numLowerBitsPerNumber;
    }

    [[nodiscard]] Index getNumUpperBitsPerNumber() const noexcept {
        return numBitsInNumber - numLowerBitsPerNumber;
    }
};

}// namespace ads

#endif//BITVECTOR_ELIAS_FANO_HPP
