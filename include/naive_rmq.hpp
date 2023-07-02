#ifndef BITVECTOR_NAIVE_RMQ_HPP
#define BITVECTOR_NAIVE_RMQ_HPP

#include "common.hpp"
#include "rand_access_iter.hpp"
#include <algorithm>

namespace ads {

#ifdef ADS_HAS_CPP20

template<typename T>
concept HasValueType = requires { typename T::value_type; };

template<typename RmqType>
struct RmqValueType : std::type_identity<Elem> {}; // The succinct rmq doesn't have a value_type

template<HasValueType RmqType>
struct RmqValueType<RmqType> : std::type_identity<typename RmqType::value_type> {};


template<typename RmqType, typename ValueType = typename RmqValueType<RmqType>::type>
concept Rmq = requires(const RmqType& r) {
    RmqType();
    RmqType(std::vector<ValueType>());
    RmqType((ValueType*)(nullptr), Index());
    RmqType(std::initializer_list<ValueType>());
    { RmqType::name } -> std::convertible_to<std::string>;
    { r(Index(), Index()) } -> std::convertible_to<Index>;
    { r.sizeInBits() } -> std::convertible_to<Index>;
};
#define ADS_RMQ_CONCEPT Rmq
#define ADS_RMQ_CONCEPT_FOR(value_type) Rmq<value_type>

#else
#define ADS_RMQ_CONCEPT class
#define ADS_RMQ_CONCEPT_FOR(value_type) class
#endif

template<typename T = Elem, typename Comparator = std::less<>>
struct [[nodiscard]] SimpleRMQ : std::vector<T> {
    using Base = std::vector<T>;

    [[no_unique_address]] Comparator comp = Comparator{};

    constexpr static const char name[] = "Simple RMQ";

    using Base::Base;

    explicit ADS_CPP20_CONSTEXPR SimpleRMQ(const std::vector<T>& vec) : Base(vec) {}

    ADS_CPP20_CONSTEXPR SimpleRMQ(T* ptr, Index length) : Base(ptr, ptr + length) {}

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rmq(Index first, Index last) const noexcept {
        assert(first < last);
        return std::min_element(this->begin() + first, this->begin() + last, comp) - this->begin();
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index operator()(Index first, Index last) const noexcept {
        return rmq(first, last);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<const T> values() const noexcept { return Span<const T>(*this); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index sizeInBits() const noexcept { return this->size() * sizeof(T) * 8; }
};


/// \brief Uses O(n^2) bits to answer range minimum queries in O(1)
/// \tparam T type of elements, ie. the `value_type`
/// \tparam Comparator Used to compare two elements, defaults to std::less (std::greater would turn this into range maximum queries)
template<typename Comparator = std::less<>>
class [[nodiscard]] NaiveRMQ {
    Allocation<Index> allocation = Allocation<Index>();
    View<Index> arr = View<Index>();
    Index length = 0;
    [[no_unique_address]] Comparator comp = Comparator{};

    // If `numValues * (numValues + 1)` overflows, we're probably in trouble anyway
    static ADS_CPP20_CONSTEXPR Index completeSize(Index numValues) noexcept {
        assert(numValues >= 0);
        return numValues * (numValues + 1) / 2;
    }

    ADS_CPP20_CONSTEXPR NaiveRMQ(Index length, CreateWithSizeTag)
        : allocation(completeSize(length)), arr(allocation.memory(), allocation.size()), length(length) {}

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index arrIndex(Index first, Index last) const noexcept {
        assert(0 <= first && first < last && last <= length);
        Index seqLen = last - 1 - first;
        Index n = length - 1 - first;
        return n * (n + 1) / 2 + seqLen;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index& minIdx(Index first, Index last) noexcept {
        return arr.ptr[arrIndex(first, last)];
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index minIdx(Index first, Index last) const noexcept {
        return arr[arrIndex(first, last)];
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index& minIdx(Index first) noexcept { return minIdx(first, first + 1); }

public:
    constexpr static const char name[] = "Naive RMQ";

    constexpr NaiveRMQ() = default;

    template<typename Range, typename = std::void_t<decltype(maybe_ranges::begin(std::declval<Range&>()))>> // no concepts in C++17
    explicit ADS_CPP20_CONSTEXPR NaiveRMQ(const Range& values)
        : NaiveRMQ(maybe_ranges::begin(values), maybe_ranges::end(values)) {}

    template<typename T>
    ADS_CPP20_CONSTEXPR NaiveRMQ(std::initializer_list<T> list) : NaiveRMQ(list.begin(), list.end()) {}

    template<typename ForwardIter, typename Sentinel>
    ADS_CPP20_CONSTEXPR NaiveRMQ(ForwardIter beginIt, Sentinel endIt) noexcept
        : NaiveRMQ(maybe_ranges::distance(beginIt, endIt), CreateWithSizeTag{}) {
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

    ADS_CPP20_CONSTEXPR NaiveRMQ& operator=(NaiveRMQ&&) noexcept = default;

    [[nodiscard]] constexpr Index rmq(Index first, Index last) const noexcept { return minIdx(first, last); }

    [[nodiscard]] constexpr Index operator()(Index first, Index last) const noexcept { return rmq(first, last); }

    [[nodiscard]] constexpr Index size() const noexcept { return length; }

    [[nodiscard]] constexpr Index sizeInBits() const noexcept { return completeSize(size()) * sizeof(Index) * 8; }

    constexpr static auto getValue = [](const auto& rmq, Index i) -> Index { return rmq.minIdx(i, i + 1); };

    using ValueIter = RandAccessIter<NaiveRMQ, decltype(getValue)>;

    [[nodiscard]] constexpr ValueIter valueIter(Index i) const { return ValueIter(*this, getValue, i); }

    [[nodiscard]] constexpr Subrange<ValueIter> values() const noexcept { return {valueIter(0), valueIter(size())}; }
};

} // namespace ads

#endif // BITVECTOR_NAIVE_RMQ_HPP
