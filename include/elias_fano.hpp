#ifndef BITVECTOR_ELIAS_FANO_HPP
#define BITVECTOR_ELIAS_FANO_HPP

#include "bitvector.hpp"
#include <cmath>
namespace ads {

template<typename Number = std::uint64_t, ADS_LAYOUT_CONCEPT BitvecLayout = SimpleLayout<>>
class EliasFano {

    Index numInts = 0;
    //    Index numLowerBitsPerNumber = 0;
    Index lowerBitMask = 0;
    Index allocatedSizeInElems = 0;
    Allocation<> allocation = Allocation<>();
    Bitvector<BitvecLayout> upper = Bitvector<BitvecLayout>();
    BitView<dynSize> lower = BitView<dynSize>();

private:
    static constexpr Index numBitsInNumber = 8 * sizeof(Number);

    EliasFano(Index numInts, CreateWithSizeTag) : numInts(numInts) {
        // add 1 before taking the log to prevent problems when numInts is 0 or 1, which would cause a shift by numBitsInNumber
        Index numLowerBitsPerNumber = numBitsInNumber - Index(std::ceil(std::log2(numInts + 1)));
        Index lowerSizeInElems = roundUpDiv(numLowerBitsPerNumber * numInts, 8 * sizeof(Elem));
        Index upperSizeInBits = numInts + (1 << (numBitsInNumber - numLowerBitsPerNumber)) + 1;
        Index allocatedUpperSizeInElems = Bitvector<BitvecLayout>::allocatedSizeInElems(roundUpDiv(upperSizeInBits, 64));
        assert(numLowerBitsPerNumber > 0 && (numLowerBitsPerNumber < numBitsInNumber || numInts == 0));

        allocatedSizeInElems = lowerSizeInElems + allocatedUpperSizeInElems;
        allocation = Allocation<>(allocatedSizeInElems);
        lowerBitMask = (Number(1) << numLowerBitsPerNumber) - 1;
        lower = BitView<dynSize>(allocation.memory(), lowerSizeInElems);
        lower.bitAccess.numBits = numLowerBitsPerNumber;
        assert(lowerSizeInElems == lower.numT);
        upper = Bitvector<BitvecLayout>(upperSizeInBits, allocation.memory() + lowerSizeInElems);
        assert(upper.allocatedSizeInElems() == allocatedUpperSizeInElems);
    }


public:
    template<typename Integer, typename = std::enable_if_t<std::is_convertible_v<Integer, Index>>>
    explicit EliasFano(Integer numInts) : EliasFano(numInts, CreateWithSizeTag{}) {}

    template<typename Range, typename = std::void_t<decltype(maybe_ranges::begin(std::declval<Range&>()))>> // no concepts in C++17
    explicit EliasFano(const Range& numbers) : EliasFano(maybe_ranges::begin(numbers), maybe_ranges::end(numbers)) {}

    EliasFano(std::initializer_list<Number> list) : EliasFano(maybe_ranges::begin(list), maybe_ranges::end(list)) {}

    // TODO: Instead of the worst case assumpption that max(first, last) == (1 << numBits) - 1, actually check the last element
    template<typename ForwardIter, typename Sentinel>
    EliasFano(ForwardIter first, Sentinel last) noexcept
        : EliasFano(maybe_ranges::distance(first, last), CreateWithSizeTag{}) {
        Elem currentUpperEntry = 0;
        Index lastUpperElemIdx = 0;
        Index elemIdx = 0;
        while (first != last) {
            Elem upperPart = Elem(*first) >> numLowerBitsPerNumber();
            Elem newBitIdx = upperPart + elemIdx + 1;
            Index currentUpperElemIdx = Index(newBitIdx / 64);
            if (currentUpperElemIdx > lastUpperElemIdx) {
                upper.setElem(lastUpperElemIdx, currentUpperEntry);
                for (Index k = lastUpperElemIdx + 1; k < currentUpperElemIdx; ++k) {
                    upper.setElem(k, 0);
                }
                lastUpperElemIdx = currentUpperElemIdx;
                currentUpperEntry = 0;
            }
            currentUpperEntry |= Elem(1) << (newBitIdx % 64); // the order of bits in an element is reversed
            lower.setBits(elemIdx, Elem(*first));
            ++elemIdx;
            ++first;
        }
        upper.setElem(lastUpperElemIdx, currentUpperEntry);
        for (++lastUpperElemIdx; lastUpperElemIdx < upper.sizeInElems(); ++lastUpperElemIdx) {
            upper.setElem(lastUpperElemIdx, 0);
        }
        upper.buildMetadata();
    }

    [[nodiscard]] const Bitvector<BitvecLayout>& getUpper() const noexcept { return upper; }

    [[nodiscard]] const BitView<dynSize>& getLower() const noexcept { return lower; }

    [[nodiscard]] Index numUpperBitsPerNumber() const noexcept { return numBitsInNumber - numLowerBitsPerNumber(); }

    [[nodiscard]] Index numLowerBitsPerNumber() const noexcept { return lower.bitAccess.numBits; }

    [[nodiscard]] Elem getUpperPart(Index i) const { return upper.selectOne(i) - i - 1; }

    [[nodiscard]] Elem getLowerPart(Index i) const { return lower.getBits(i); }

    [[nodiscard]] Number get(Index i) const {
        if (i < 0 || i >= size()) [[unlikely]] {
            throw std::invalid_argument("EliasFano::get() Index out of range");
        }
        Elem upperPart = getUpperPart(i);
        Elem lowerPart = getLowerPart(i);
        return Number(lowerPart + (upperPart << numLowerBitsPerNumber()));
    }

    [[nodiscard]] Number operator[](Index i) const { return get(i); }

    constexpr static auto getNumber = [](const auto& ef, Index i) -> Number { return ef.get(i); };

    using NumberIter = RandAccessIter<EliasFano, decltype(getNumber)>;

    NumberIter numberIter(Index i) const { return NumberIter(*this, getNumber, i); }

    Subrange<NumberIter> numbers() const noexcept { return {numberIter(0), numberIter(size())}; }

    [[nodiscard]] Number predecessor(Number n) const {
        Number upperSearchBits = n >> numLowerBitsPerNumber();
        Number lowerSearchBits = n & lowerBitMask;
        Index first = upper.selectZero(upperSearchBits) - upperSearchBits;
        Index last = upper.selectZero(upperSearchBits + 1) - upperSearchBits - 1;
        for (Index i = last - 1; i >= first; --i) {
            Number lowerBits = lower.getBits(i);
            if (lowerBits <= lowerSearchBits) {
                return lowerBits + (upperSearchBits << numLowerBitsPerNumber());
            }
        }
        // there was no element with the same upper part, so return the greatest element with a smaller upper part
        if (first == 0) [[unlikely]] {
            throw std::invalid_argument("no predecessor found for " + std::to_string(n));
        }
        return get(first - 1);
    }

    [[nodiscard]] Number successor(Number n) const {
        Number upperSearchBits = n >> numLowerBitsPerNumber();
        Number lowerSearchBits = n & lowerBitMask;
        Index first = upper.selectZero(upperSearchBits) - upperSearchBits;
        Index last = upper.selectZero(upperSearchBits + 1) - upperSearchBits - 1;
        for (Index i = first; i != last; ++i) {
            Number lowerBits = lower.getBits(i);
            if (lowerBits >= lowerSearchBits) {
                return lowerBits + (upperSearchBits << numLowerBitsPerNumber());
            }
        }
        // there was no element with the same upper part, so return the greatest element with a smaller upper part
        if (last == numInts) [[unlikely]] {
            throw std::invalid_argument("no successor found for " + std::to_string(n));
        }
        return get(last);
    }

    [[nodiscard]] Index size() const noexcept { return numInts; }

    [[nodiscard]] Index numAllocatedBits() const noexcept { return allocatedSizeInElems * sizeof(Elem) * 8; }
};

} // namespace ads

#endif // BITVECTOR_ELIAS_FANO_HPP
