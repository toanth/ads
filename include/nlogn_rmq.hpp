#ifndef BITVECTOR_NLOGN_RMQ_HPP
#define BITVECTOR_NLOGN_RMQ_HPP

#include "bit.hpp"
#include "bit_access.hpp"
#include "common.hpp"
namespace ads {

ADS_CPP20_CONSTEXPR Index numLevels(Index length) noexcept {
    assert(length >= 0);
    if (length <= 1) {
        return length;
    }
    return log2(Elem(length)) + 1;
}

ADS_CPP20_CONSTEXPR Index minimaSize(Index length) noexcept {
    Index levels = numLevels(length);
    return length * levels; // TODO: Save some elements here?
}


/// \brief \Operations for the  n log n rmq structure. A separate class template because the linear rmq also uses a
/// slightly different n log n rmq, which is also implemented using this class template.
/// \tparam Derived The actual RMQ type.
/// \tparam Comp comparator, defaults to std::less<>
template<typename Derived, typename Comp>
class NLogNRmqOps {

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

    Index operator()(Index lower, Index upper) const { return rmq(lower, upper); }

    [[nodiscard]] Index size() const noexcept { return derived().length; }

    [[nodiscard]] Comp getComp() const noexcept { return comp; }
};


template<typename T = Elem, typename IndexType = Index, typename Comp = std::less<>>
struct NLogNRmq : NLogNRmqOps<NLogNRmq<T, IndexType, Comp>, Comp> {
    using Base = NLogNRmqOps<NLogNRmq<T, IndexType, Comp>, Comp>;
    friend Base;
    Allocation<T> allocation = Allocation<T>(); // must be the first data member so that its destructor is called last
    View<T> array = View<T>();
    View<IndexType> minima = View<IndexType>();
    Index length = 0;

    static_assert(alignof(T) % alignof(IndexType) == 0);

    constexpr static const char name[] = "n log n space RMQ (1)";

    NLogNRmq() = default;

    NLogNRmq(std::initializer_list<T> list) : NLogNRmq(Span<const T>(list.begin(), list.end())) {}

    NLogNRmq(Index length, CreateWithSizeTag, T* mem = nullptr)
        : allocation(completeSize(length), mem), length(length) {
        T* ptr = allocation.memory();
        array = View<T>(ptr, length);
        minima = View<IndexType>(ptr + length, completeSize(length) - length);
    }

    explicit NLogNRmq(Span<const T> values, T* mem = nullptr) : NLogNRmq(values.size(), CreateWithSizeTag{}, mem) {
        this->init(values);
    }

    explicit NLogNRmq(const std::vector<T>& values, T* mem = nullptr) : NLogNRmq(Span<const T>(values), mem) {}

    NLogNRmq(T* ptr, Index length, T* mem = nullptr) : NLogNRmq(Span<const T>(ptr, length), mem) {}


    static Index completeSize(Index length) noexcept {
        return length + roundUpDiv(minimaSize(length) * sizeof(IndexType), sizeof(T));
    }

    [[nodiscard]] Index sizeInBits() const noexcept { return completeSize(length) * sizeof(T) * 8; }

    const T& operator[](Index i) const noexcept { return getArrayElement(i); }

    [[nodiscard]] Span<const T> values() const noexcept { return Span<const T>(array.ptr.get(), this->size()); }

private:
    T& getArrayElement(Index i) noexcept { return array.ptr[i]; }
    const T& getArrayElement(Index i) const noexcept { return array.ptr[i]; }


    IndexType getMinimum(Index i) const noexcept { return minima[i]; }
    void setMinimum(Index i, const IndexType& value) noexcept { minima.setBits(i, value); }
};

} // namespace ads

#endif // BITVECTOR_NLOGN_RMQ_HPP
