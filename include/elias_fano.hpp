#ifndef BITVECTOR_ELIAS_FANO_HPP
#define BITVECTOR_ELIAS_FANO_HPP

#include "bitvector/recursive_bitvec.hpp" // TODO: Use better bitvector implementation
#include <cmath>
namespace ads {

/// \brief Elias-Fano which also supports predecessor and successor queries. The upper bitvector represents  occurring
/// values with zeros instead of ones, ie it inverts all bits, which allows using selectOne() instead of selectZero() in
/// the common case.
/// \tparam Number The type of numbers stored the Elias-Fano data structure. Currently, this must be a built-in integer.
/// \tparam Bitvec The bitvector to use for the upper part.
template<typename Number = std::uint64_t, ADS_BITVEC_CONCEPT Bitvec = EfficientSelectBitvec<>>
class [[nodiscard]] EliasFano {
    static_assert(std::is_integral_v<Number>, "Elias fano only works for integers");

    Index numInts = 0;
    Index lowerBitMask = 0;
    Index allocatedSizeInBytes = 0;
    Elem smallestNumber = 0;
    Elem largestNumber = 0;
    Index bitsPerNumber = sizeof(Number) * 8; // usually overwritten in the ctor
    Allocation<> allocation = Allocation<>();
    Bitvec upper = Bitvec();
    BitView<dynSize> lower = BitView<dynSize>();
    constexpr static Index linearFallbackSize = 8;

    template<typename Range>
    ADS_CPP20_CONSTEXPR void build(const Range& range) noexcept {
        Elem currentUpperEntry = 0;
        Index lastUpperElemIdx = 0;
        Index i = 0;
        for (Number valUnmodified : range) {
            Elem val = Elem(valUnmodified) - smallestNumber;
            //            ADS_ASSUME(valUnmodified >= smallestNumber);
            Index upperPart = Index(val >> numLowerBitsPerNumber());
            Index newBitIdx = upperPart + i + 1;
            Index currentUpperElemIdx = newBitIdx / 64;
            if (currentUpperElemIdx > lastUpperElemIdx) {
                upper.setLimb(lastUpperElemIdx, ~currentUpperEntry);
                for (Index k = lastUpperElemIdx + 1; k < currentUpperElemIdx; ++k) {
                    upper.setLimb(k, Elem(-1));
                }
                lastUpperElemIdx = currentUpperElemIdx;
                currentUpperEntry = 0;
            }
            currentUpperEntry |= Elem(1) << (newBitIdx % 64); // the order of bits in an element is reversed
            if (numLowerBitsPerNumber() > 0) {
                lower.setBits(i, Elem(val));
            }
            ++i;
        }
        upper.setLimb(lastUpperElemIdx, ~currentUpperEntry);
        for (++lastUpperElemIdx; lastUpperElemIdx < upper.numLimbs(); ++lastUpperElemIdx) {
            upper.setLimb(lastUpperElemIdx, Elem(-1));
        }
        upper.buildMetadata();
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Elem predecessorImpl(Elem n) const {
        Number upperSearchBits = n >> numLowerBitsPerNumber();
        Number lowerSearchBits = n & lowerBitMask;
        ADS_ASSUME(upper.numOnes() > upperSearchBits);
        auto [first, last] = upper.selectOneAndPrevOne(upperSearchBits + 1);
        first -= upperSearchBits;
        last -= upperSearchBits + 1;
        if (numLowerBitsPerNumber() == 0) {
            return first == last ? getImpl(first - 1) : getImpl(first);
        }
        while (last - first > linearFallbackSize) {
            Index mid = (first + last) / 2;
            Number lowerBits = lower.getBits(mid);
            if (lowerBits > lowerSearchBits) {
                last = mid;
            } else {
                first = mid;
            }
        }
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
        return getImpl(first - 1);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Elem successorImpl(Elem n) const {
        Number upperSearchBits = n >> numLowerBitsPerNumber();
        Number lowerSearchBits = n & lowerBitMask;
        auto [first, last] = upper.selectOneAndPrevOne(upperSearchBits + 1);
        first -= upperSearchBits;
        last -= upperSearchBits + 2;
        if (numLowerBitsPerNumber() == 0) {
            return getImpl(first);
        }
        while (last - first > linearFallbackSize) {
            Index mid = (last + first) / 2;
            Number lowerBits = lower.getBits(mid);
            if (lowerBits <= lowerSearchBits) {
                first = mid;
            } else {
                last = mid;
            }
        }
        for (Index i = first; i <= last; ++i) {
            Number lowerBits = lower.getBits(i);
            if (lowerBits >= lowerSearchBits) {
                return lowerBits + (upperSearchBits << numLowerBitsPerNumber());
            }
        }
        // there was no element with the same upper part, so return the greatest element with a smaller upper part
        if (last == numInts - 1) [[unlikely]] {
            throw std::invalid_argument("no successor found for " + std::to_string(n));
        }
        return getImpl(last + 1);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Elem getImpl(Index i) const {
        if (i < 0 || i >= size()) {
            throw std::invalid_argument("EliasFano::get() Index out of range");
        }
        Elem upperPart = getUpperPart(i);
        if (numLowerBitsPerNumber() == 0) {
            return upperPart;
        }
        Elem lowerPart = getLowerPart(i);
        return lowerPart + (upperPart << numLowerBitsPerNumber());
    }

    /// \brief Use the bitvector to find the upper parts of the ith stored number.
    [[nodiscard]] constexpr Elem getUpperPart(Index i) const { return upper.selectZero(i) - i - 1; }

    /// \brief The ith entry in the lower array, which holds the lower bits of the ith stored number.
    [[nodiscard]] constexpr Elem getLowerPart(Index i) const { return lower.getBits(i); }

public:
    template<typename Range, typename = std::void_t<decltype(maybe_ranges::begin(std::declval<Range&>()))>> // no concepts in C++17
    explicit ADS_CPP20_CONSTEXPR EliasFano(const Range& numbers) {
        numInts = maybe_ranges::size(numbers);
        Elem rangeOfValues = 0;
        if (numInts > 0) [[likely]] {
            auto iter = maybe_ranges::begin(numbers);
            smallestNumber = Elem(*iter);
            maybe_ranges::advance(iter, numInts - 1);
            largestNumber = *iter;
            rangeOfValues = largestNumber - smallestNumber;
        }

        bitsPerNumber = 1 + intLog2(std::max(rangeOfValues, Elem(1)));
        Index upperBitsPerNumber = roundUpLog2(Elem(numInts + 2)); // + 2 to prevent UB for numInts <= 1
        bitsPerNumber = std::max(bitsPerNumber, upperBitsPerNumber);
        Index lowerBitsPerNumber = bitsPerNumber - upperBitsPerNumber;
        Index lowerSizeInBytes = roundUpTo(lowerBitsPerNumber * numInts, 8 * CACHELINE_SIZE_BYTES) / 8;
        Index maxUpperVal = Index(rangeOfValues >> lowerBitsPerNumber);
        Index numBitsInUpper = 2 + numInts + maxUpperVal; // +2 because the first bit is always a one, and the last bit as well
        Index allocatedUpperSizeInBytes = Bitvec::allocatedSizeInBytesForBits(numBitsInUpper);
        ADS_ASSUME(allocatedUpperSizeInBytes % CACHELINE_SIZE_BYTES == 0);
        ADS_ASSUME(maxUpperVal < Index(1) << upperBitsPerNumber);

        allocatedSizeInBytes = lowerSizeInBytes + allocatedUpperSizeInBytes;
        allocation = Allocation<>(allocatedSizeInBytes);
        auto* memory = allocation.memory();
        ADS_ASSUME_ALIGNED(memory, 32);
        upper = Bitvec::uninitializedForSize(numBitsInUpper, memory);
        ADS_ASSUME(upper.allocatedSizeInBytes() == allocatedUpperSizeInBytes);
        lower = BitView<dynSize>(upper.alloc().endOfMemory(), lowerSizeInBytes / sizeof(Limb));
        lowerBitMask = (Number(1) << lowerBitsPerNumber) - 1;
        lower.bitAccess.numBits = lowerBitsPerNumber;
        ADS_ASSUME(lowerSizeInBytes == lower.numT * sizeof(Limb));
        build(numbers);
    }

    ADS_CPP20_CONSTEXPR EliasFano(std::initializer_list<Number> list) noexcept
        : EliasFano(maybe_ranges::begin(list), maybe_ranges::end(list)) {}

    template<typename ForwardIter, typename Sentinel>
    ADS_CPP20_CONSTEXPR EliasFano(ForwardIter first, Sentinel last) noexcept
        : EliasFano(Subrange<ForwardIter, Sentinel>{first, last}) {}

    [[nodiscard]] constexpr Index numBitsPerNumber() const noexcept { return bitsPerNumber; }

    [[nodiscard]] constexpr Number get(Index i) const { return Number(getImpl(i) + smallestNumber); }

    [[nodiscard]] constexpr Number operator[](Index i) const { return get(i); }

    constexpr static auto getNumber = [](const auto& ef, Index i) -> Number { return ef.get(i); };

    using NumberIter = RandAccessIter<EliasFano, decltype(getNumber)>;

    constexpr NumberIter numberIter(Index i) const { return NumberIter(*this, getNumber, i); }

    constexpr Subrange<NumberIter> numbers() const noexcept { return {numberIter(0), numberIter(size())}; }

    [[nodiscard]] [[using gnu: hot, pure]] constexpr Number predecessor(Number n) const {
        // the cast is implementation defined before C++20 (but still does the right thing)
        if (n >= I64(largestNumber)) {
            return Number(largestNumber);
        }
        return Number(predecessorImpl(Elem(n) - smallestNumber) + smallestNumber);
    }


    [[nodiscard]] constexpr Number successor(Number n) const {
        // the cast is implementation defined before C++20 (but still does the right thing)
        if (n <= I64(smallestNumber)) {
            return Number(smallestNumber);
        }
        return Number(successorImpl(Elem(n) - smallestNumber) + smallestNumber);
    }

    /// \brief Returns the smallest value, which is the same as `*numbers().begin()` but more efficient
    [[nodiscard]] constexpr Number getSmallest() const noexcept { return Number(smallestNumber); }

    /// \brief Returns the largest values, which is the same as *(numbers.end() - 1) but more efficient
    [[nodiscard]] constexpr Number getLargest() const noexcept { return Number(largestNumber); }

    /// \brief The number of stored values.
    [[nodiscard]] constexpr Index size() const noexcept { return numInts; }

    /// \brief The total amount of bits allocated on the heap by this class, including its data members.
    /// Note that data members within this class itself, such as pointers to allocated memory, don't count.
    [[nodiscard]] constexpr Index numAllocatedBits() const noexcept {
        ADS_ASSUME(allocatedSizeInBytes == upper.allocatedSizeInBytes() + lower.sizeInBytes());
        return allocatedSizeInBytes * 8;
    }

    [[nodiscard]] constexpr const Bitvec& getUpper() const noexcept { return upper; }

    [[nodiscard]] constexpr const BitView<dynSize>& getLower() const noexcept { return lower; }

    [[nodiscard]] constexpr Index numUpperBitsPerNumber() const noexcept {
        return numBitsPerNumber() - numLowerBitsPerNumber();
    }

    [[nodiscard]] constexpr Index numLowerBitsPerNumber() const noexcept { return lower.bitAccess.numBits; }
};


} // namespace ads

#endif // BITVECTOR_ELIAS_FANO_HPP
