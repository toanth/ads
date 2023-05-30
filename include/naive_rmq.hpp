#ifndef BITVECTOR_NAIVE_RMQ_HPP
#define BITVECTOR_NAIVE_RMQ_HPP

#include "common.hpp"
#include "rand_access_iter.hpp"
#include <algorithm>

namespace ads {

template<typename T = Elem, typename Comparator = std::less<>>
struct SimpleRMQ : std::vector<T> {
    using Base = std::vector<T>;

    [[no_unique_address]] Comparator comp = Comparator{};

    using Base::Base;

    SimpleRMQ(const std::vector<T>& vec) : Base(vec) {}

    [[nodiscard]] Index rmq(Index first, Index last) const noexcept {
        assert(first < last);
        return std::min_element(this->begin() + first, this->begin() + last, comp) - this->begin();
    }

    [[nodiscard]] Index operator()(Index first, Index last) const noexcept {
        return rmq(first, last);
    }
};


/// \brief Uses O(n^2) bits to answer range minimum queries in O(1)
/// \tparam T type of elements, ie. the `value_type`
/// \tparam Comparator Used to compare two elements, defaults to std::less (std::greater would turn this into range maximum queries)
template<typename Comparator = std::less<>>
class NaiveRMQ {
    std::unique_ptr<Index[]> arr;
    Index length;
    [[no_unique_address]] Comparator comp = Comparator{};


    // If `length * (length + 1)` overflows, we're probably in trouble anyway
    NaiveRMQ(Index length, CreateWithSizeTag) : arr(makeUniqueForOverwrite<Index>(length * (length + 1) / 2)), length(length) {
    }

    Index arrIndex(Index first, Index last) const noexcept {
        assert(0 <= first && first < last && last <= length);
        Index seqLen = last - 1 - first;
        Index n = length - 1 - first;
        return n * (n + 1) / 2 + seqLen;
    }

    Index& minIdx(Index first, Index last) noexcept {
        return arr[arrIndex(first, last)];
    }
    Index minIdx(Index first, Index last) const noexcept {
        return arr[arrIndex(first, last)];
    }
    Index& minIdx(Index first) noexcept {
        return minIdx(first, first + 1);
    }

public:
    template<typename Range, typename = std::void_t<decltype(maybe_ranges::begin(std::declval<Range&>()))>>// no concepts in C++17
    explicit NaiveRMQ(const Range& values)
        : NaiveRMQ(maybe_ranges::begin(values), maybe_ranges::end(values)) {
    }

    template<typename T>
    NaiveRMQ(std::initializer_list<T> list) : NaiveRMQ(list.begin(), list.end()) {}

    template<typename ForwardIter, typename Sentinel>
    NaiveRMQ(ForwardIter beginIt, Sentinel endIt) noexcept : NaiveRMQ(maybe_ranges::distance(beginIt, endIt), CreateWithSizeTag{}) {
        Index i = 0;
        std::vector<Elem> values(length);
        for (auto iter = beginIt; iter != endIt; ++iter, ++i) {
            values[i] = *iter;
            minIdx(i) = i;
        }
        for (Index first = 0; first < length; ++first) {
            for (Index last = first + 1; last < length; ++last) {
                Index m = minIdx(first, last);
                minIdx(first, last + 1) = comp(values[m], values[last]) ? m : last;
            }
        }
    }

    [[nodiscard]] Index rmq(Index first, Index last) const noexcept {
        return minIdx(first, last);
    }

    [[nodiscard]] Index operator()(Index first, Index last) const noexcept {
        return rmq(first, last);
    }

    [[nodiscard]] Index size() const noexcept {
        return length;
    }

    constexpr static auto getValue = [](const NaiveRMQ& rmq, Index i) -> Index {
        return rmq.minIdx(i, i + 1);
    };

    using ValueIter = RandAccessIter<NaiveRMQ, decltype(getValue)>;

    [[nodiscard]] ValueIter valueIter(Index i) const {
        return ValueIter(*this, getValue, i);
    }

    [[nodiscard]] Subrange<ValueIter> values() const noexcept {
        return {valueIter(0), valueIter(size())};
    }
};

}// namespace ads

#endif//BITVECTOR_NAIVE_RMQ_HPP
