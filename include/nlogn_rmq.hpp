#ifndef BITVECTOR_NLOGN_RMQ_HPP
#define BITVECTOR_NLOGN_RMQ_HPP

#include "common.hpp"
namespace ads {

ADS_CPP20_CONSTEXPR Index numLevels(Index length) noexcept {
    assert(length >= 0);
    if (length <= 1) { return length; }
    return log2(Elem(length)) + 1;
}

ADS_CPP20_CONSTEXPR Index minimaSize(Index length) noexcept {
    Index levels = numLevels(length);
    return length * levels;// TODO: Save some elements here?
}


template<typename Derived, typename Comp>
class NlognRmqOps {

    [[no_unique_address]] Comp comp = Comp();

    const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    Derived& derived() noexcept { return static_cast<Derived&>(*this); }

    //    auto& getArray() noexcept { return derived().array; }
    //    const auto& getArray() const noexcept { return derived().array; }

public:
    void init() noexcept {
        for (Index i = 0; i < size(); ++i) {
            derived().setMinimum(minimaIndex(i, 0), i);
        }
        Index levels = numLevels(size());
        for (Index lvl = 1; lvl < levels; ++lvl) {
            for (Index i = 0, j = Index(1) << lvl; j <= size(); ++i, ++j) {
                Index leftMinIdx = getMinimumStartingAt(i, lvl - 1);
                Index rightMinIdx = getMinimumEndingAt(j, lvl - 1);
                Index minIdx = comp(derived().getArrayElement(leftMinIdx), derived().getArrayElement(rightMinIdx)) ? leftMinIdx : rightMinIdx;
                derived().setMinimum(minimaIndex(i, lvl), minIdx);
            }
        }
    }

    template<typename Container>
    void init(const Container& values, Comp c = Comp()) noexcept {
        comp = c;
        for (Index i = 0; i < values.size(); ++i) {
            derived().getArrayElement(i) = values[i];
        }
        init();
    }

    [[nodiscard]] Index minimaIndex(Index start, Index log2Length) const noexcept {
        assert(start >= 0 && start < size() && log2Length >= 0 && log2Length < numLevels(size()));
        return log2Length * size() + start;
    }

    [[nodiscard]] Index getMinimumStartingAt(Index start, Index log2Length) const noexcept {
        return derived().getMinimum(minimaIndex(start, log2Length));
    }

    [[nodiscard]] Index getMinimumEndingAt(Index end, Index log2Length) const noexcept {
        return getMinimumStartingAt(end - (Index(1) << log2Length), log2Length);
    }

    [[nodiscard]] Index rmq(Index lower, Index upper) const {
        if (!(0 <= lower && lower < upper && upper <= size())) [[unlikely]] {
            throw std::invalid_argument("Invalid rmq range");
        }
        //        if (lower == upper) { return lower; }
        Index log2Length = log2(Elem(upper - lower));
        Index leftMin = getMinimumStartingAt(lower, log2Length);
        Index rightMin = getMinimumEndingAt(upper, log2Length);
        return comp(derived().getArrayElement(leftMin), derived().getArrayElement(rightMin)) ? leftMin : rightMin;
    }

    Index operator()(Index lower, Index upper) const {
        return rmq(lower, upper);
    }

    [[nodiscard]] Index size() const noexcept { return derived().length; }

    [[nodiscard]] Comp getComp() const noexcept { return comp; }
};


template<typename T, typename IndexType, typename Comp>
struct NlognRMQSeparateArrays : NlognRmqOps<NlognRMQSeparateArrays<T, IndexType, Comp>, Comp> {
    using Base = NlognRmqOps<NlognRMQSeparateArrays<T, IndexType, Comp>, Comp>;
    friend Base;
    std::unique_ptr<T[]> array = nullptr;
    std::unique_ptr<IndexType[]> minima = nullptr;
    Index length = 0;

    constexpr static const char name[] = "n log n space RMQ (2)";

    NlognRMQSeparateArrays() = default;

    NlognRMQSeparateArrays(std::initializer_list<T> list) : NlognRMQSeparateArrays(Span<const T>(list)) {}

    NlognRMQSeparateArrays(Index length, CreateWithSizeTag)
        : array(makeUniqueForOverwrite<T>(length)), minima(makeUniqueForOverwrite<IndexType>(minimaSize(length))), length(length) {}

    explicit NlognRMQSeparateArrays(Span<const T> values) : NlognRMQSeparateArrays(values.size(), CreateWithSizeTag{}) {
        this->init(values);
    }

    const T& operator[](Index i) const noexcept { return getArrayElement(i); }

private:
    T& getArrayElement(Index i) noexcept { return array[i]; }
    const T& getArrayElement(Index i) const noexcept { return array[i]; }

    [[nodiscard]] Index getMinimum(Index i) const noexcept {
        return minima[i];
    }
    void setMinimum(Index i, const IndexType& value) const noexcept {
        minima[i] = value;
    }
};

template<typename T, typename IndexType, typename Comp>
struct NlognRMQOneArray : NlognRmqOps<NlognRMQOneArray<T, IndexType, Comp>, Comp> {
    using Base = NlognRmqOps<NlognRMQOneArray<T, IndexType, Comp>, Comp>;
    friend Base;
    std::unique_ptr<T[]> array = nullptr;
    unsigned char* ADS_RESTRICT minima = nullptr;
    Index length = 0;

    constexpr static const char name[] = "n log n space RMQ (1)";

    NlognRMQOneArray() = default;

    NlognRMQOneArray(std::initializer_list<T> list) : NlognRMQOneArray(Span<const T>(list.begin(), list.end())) {}

    NlognRMQOneArray(Index length, CreateWithSizeTag) : array(makeUniqueForOverwrite<T>(completeSize(length))), length(length) {
        minima = reinterpret_cast<unsigned char*>(array.get() + length);
    }

    explicit NlognRMQOneArray(Span<const T> values) : NlognRMQOneArray(values.size(), CreateWithSizeTag{}) {
        this->init(values);
    }

    static Index completeSize(Index length) noexcept {
        return length + minimaSize(length);
    }

    const T& operator[](Index i) const noexcept { return getArrayElement(i); }

private:
    T& getArrayElement(Index i) noexcept { return array[i]; }
    const T& getArrayElement(Index i) const noexcept { return array[i]; }


    IndexType getMinimum(Index i) const noexcept {
        return ptrBitCast<Index>(minima + i * sizeof(Index));
    }
    void setMinimum(Index i, const IndexType& value) const noexcept {
        std::memcpy(minima + i * sizeof(T), &value, sizeof(T));
    }
};

template<typename T = Elem, typename IndexType = Index, typename Comp = std::less<>>
using NlognRMQ = std::conditional_t<sizeof(T) == sizeof(IndexType), NlognRMQOneArray<T, IndexType, Comp>, NlognRMQSeparateArrays<T, IndexType, Comp>>;


}// namespace ads

#endif//BITVECTOR_NLOGN_RMQ_HPP
