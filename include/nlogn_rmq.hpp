#ifndef BITVECTOR_NLOGN_RMQ_HPP
#define BITVECTOR_NLOGN_RMQ_HPP

#include "common.hpp"
namespace ads {

ADS_CPP20_CONSTEXPR Index numLevels(Index length) noexcept {
    return log2(Elem(length - 1)) + 1;
}

ADS_CPP20_CONSTEXPR Index minimaSize(Index length) noexcept {
    if (length <= 1) { return 0; }
    Index levels = numLevels(length);
    return length * levels;// TODO: Save some elements here?
}

template<typename T, typename IndexType, bool>
struct NlognRMQBase {
    std::unique_ptr<T[]> array;
    std::unique_ptr<IndexType[]> minima;
    Index length;

    NlognRMQBase(Index length, CreateWithSizeTag) : array(makeUniqueForOverwrite<T>(length)), minima(minimaSize(length)), length(length) {}

    Index getMinimum(Index i) const noexcept {
        return minima[i];
    }
    void setMinimum(Index i, const IndexType& value) const noexcept {
        minima[i] = value;
    }
};

template<typename T, typename IndexType>
struct NlognRMQBase<T, IndexType, true> {
    std::unique_ptr<T[]> array;
    unsigned char* ADS_RESTRICT minima = nullptr;
    Index length;

    NlognRMQBase(Index length, CreateWithSizeTag) : array(makeUniqueForOverwrite<T>(completeSize(length))), length(length) {
        minima = reinterpret_cast<unsigned char*>(array.get() + length);
    }


    static Index completeSize(Index length) noexcept {
        return length + minimaSize(length);
    }

    IndexType getMinimum(Index i) const noexcept {
        return ptrBitCast<Index>(minima + i * sizeof(Index));
    }
    void setMinimum(Index i, const IndexType& value) const noexcept {
        std::memcpy(minima + i * sizeof(T), &value, sizeof(T));
    }
};

// Better idea(?): Don't store original array, only store results of comparison against neighbouring minima of same size? Dosn't work?
template<typename T = Elem, typename IndexType = Index, typename Comparator = std::less<>>
class NlognRMQ : NlognRMQBase<T, IndexType, sizeof(T) == sizeof(IndexType)> {

    using Base = NlognRMQBase<T, IndexType, sizeof(T) == sizeof(IndexType)>;
    using Base::array;
    using Base::length;
    using Base::minima;
    [[no_unique_address]] Comparator comp = Comparator();

public:
    using Base::getMinimum;
    using Base::setMinimum;

    NlognRMQ(Span<const T> values, Comparator comp = Comparator()) noexcept : Base(values.size(), CreateWithSizeTag{}), comp(comp) {
        std::copy(values.begin(), values.end(), array.get());
        for (Index i = 0; i < size(); ++i) {
            setMinimum(minimaIndex(i, 0), i);
        }
        Index levels = numLevels(size());
        for (Index lvl = 1; lvl < levels; ++lvl) {
            for (Index i = 0, j = IndexType(1) << lvl; j <= size(); ++i, ++j) {
                IndexType leftMinIdx = getMinimum(i, lvl - 1);
                IndexType rightMinIdx = getMinimumEndingAt(j, lvl - 1);
                IndexType minIdx = comp(values[leftMinIdx], values[rightMinIdx]) ? leftMinIdx : rightMinIdx;
                setMinimum(minimaIndex(i, lvl), minIdx);
            }
        }
    }

    [[nodiscard]] Index minimaIndex(Index start, Index log2Length) const noexcept {
        assert(start >= 0 && start < size() && log2Length >= 0 && log2Length < numLevels(size()));
        return log2Length * size() + start;
    }

    IndexType getMinimum(IndexType start, IndexType log2Length) const noexcept {
        return getMinimum(minimaIndex(start, log2Length));
    }

    IndexType getMinimumEndingAt(IndexType end, Index log2Length) const noexcept {
        return getMinimum(end - (IndexType(1) << log2Length), log2Length);
    }

    Index rmq(IndexType lower, IndexType upper) const {
        if (!(0 <= lower && lower < upper && upper <= size())) [[unlikely]] {
            throw std::invalid_argument("Invalid rmq range");
        }
        if (lower == upper) { return lower; }
        IndexType log2Size = log2(Elem(upper - lower));
        IndexType leftMin = getMinimum(lower, log2Size);
        IndexType rightMin = getMinimumEndingAt(upper, log2Size);
        return array[leftMin] < array[rightMin] ? leftMin : rightMin;
    }

    Index operator()(IndexType lower, IndexType upper) const {
        return rmq(lower, upper);
    }

    [[nodiscard]] Index size() const noexcept { return length; }
};


}// namespace ads

#endif//BITVECTOR_NLOGN_RMQ_HPP
