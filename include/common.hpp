#ifndef ADS_COMMON_HPP
#define ADS_COMMON_HPP

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>

static_assert(__cplusplus >= 201703L, "This library requires at least C++17 (some features require C++20)");
static_assert(sizeof(void*) >= 8, "This library requires a 64 bit system");

#ifdef _MSC_VER

#define ADS_HAS_MSVC_INTRINSICS
#include <intrin.h>

#elif defined __clang__ || defined __GNUC__

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++20-attribute-extensions" // else, clang complains in C++17 mode about [[likely]]
// else, clang complains about function calls in ADS_ASSUME, but they are useful in debug mode when ADS_ASSUME is assert
#pragma clang diagnostic ignored "-Wassume"
#endif

#define ADS_HAS_DEFAULT_GCC_INTRINSICS

#if defined __x86_64__
#define ADS_HAS_GCC_X86_INTRINSICS
#include "x86intrin.h"
#ifdef __BMI2__
#define ADS_HAS_GCC_BMI2
#endif
#endif

#endif // _MSC_VER


#if __cplusplus >= 202002L
#include <bit>
#include <ranges>
#endif

namespace ads {

#if __cplusplus >= 202002L
constexpr static bool hasCpp20 = true;
#define ADS_HAS_CPP20
#define ADS_CONSTEVAL consteval
#define ADS_CPP20_CONSTEXPR constexpr
#define ADS_INTEGRAL std::integral

#if __cplusplus >= 202300L
#define ADS_IF_CONSTEVAL if consteval
#else
#define ADS_IF_CONSTEVAL if (std::is_constant_evaluated())
#endif

namespace maybe_ranges = std::ranges;

#else

constexpr static bool hasCpp20 = false;
#define ADS_CONSTEVAL constexpr
#define ADS_CPP20_CONSTEXPR inline
#define ADS_INTEGRAL typename

namespace maybe_ranges = std;

#define ADS_IF_CONSTEVAL if constexpr (false)

#endif

// TODO: Make unsigned? Makes division, modulo by powers of two more efficient
using Index = std::ptrdiff_t;
using U64 = std::uint64_t;
using I64 = std::int64_t;
// TODO: Instead of a global Elem alias, define per template to allow smaller sizes
using Elem = U64;



#if __has_cpp_attribute(assume)
#define ADS_ASSUME_IMPL(x) [[assume(x)]]
#elif defined __clang__
#define ADS_ASSUME_IMPL(x) __builtin_assume(x)
#elif defined(__GNUC__) && !defined(__ICC)
// the following doesn't ignore side effects, but unfortunately, there is no better way to implement this in gcc
// -- since assume is only used internally, this isn't a huge deal
#define ADS_ASSUME_IMPL(x)                                                                                             \
    if (x) {                                                                                                           \
    } else {                                                                                                           \
        __builtin_unreachable();                                                                                       \
    }
#elif defined _MSC_VER || defined __ICC
#define ADS_ASSUME_IMPL(x) __assume(x)
#endif // __has_cpp_attribute(assume)

#ifdef NDEBUG
#define ADS_ASSUME(x)                                                                                                  \
    assert(x); /* will probably be optimized to assert(true), but still useful on the off-chance it won't */           \
    ADS_ASSUME_IMPL(x)
#else // don't use assume as that could (and probably would) cause the compiler to treat the assert as assert(true)
#define ADS_ASSUME(x) assert(x)
#endif

struct CreateWithSizeTag {};

namespace detail {

template<typename T>
struct TypeIdentity { // std::type_identity is a C++20 feature
    using Type = T;
};

template<Index NumBytes>
struct IntTypeImpl : TypeIdentity<std::uint64_t> {
    static_assert(NumBytes > 4 && NumBytes <= 8);
};

template<>
struct IntTypeImpl<4> : TypeIdentity<std::uint32_t> {};
template<>
struct IntTypeImpl<3> : TypeIdentity<std::uint32_t> {};
template<>
struct IntTypeImpl<2> : TypeIdentity<std::uint16_t> {};
template<>
struct IntTypeImpl<1> : TypeIdentity<std::uint8_t> {};

} // namespace detail

template<Index NumBytes>
using IntType = typename detail::IntTypeImpl<NumBytes>::Type;

#define ADS_RESTRICT __restrict


template<typename T, typename... Args>
ADS_CPP20_CONSTEXPR T* constructAt(T* location, Args&&... args) {
#ifdef ADS_HAS_CPP20
    return std::construct_at(location, std::forward<Args>(args)...);
#else
    return ::new (static_cast<void*>(location)) T(std::forward<Args>(args)...);
#endif
}

template<typename T>
ADS_CPP20_CONSTEXPR void uninitializedValueConstructN(T* ADS_RESTRICT dest, Index n) noexcept {
    assert(n >= 0);
    ADS_IF_CONSTEVAL {
        for (Index i = 0; i < n; ++i) {
            constructAt(dest + i);
        }
    }
    else {
        std::uninitialized_value_construct_n(dest, n);
    }
}

/// Unlike std::from_chars, the compile time version doesn't return an error code on errors,
/// instead throwing an exception which results in a compile time error.
/// Also, this version only handles non-negative integers with base <= 16 and requires the entire input to be valid.
template<ADS_INTEGRAL T>
constexpr std::from_chars_result fromChars(const char* first, const char* last, T& value, int base = 10) {
    ADS_IF_CONSTEVAL {
        value = T(0);
        T digit;
        const char* it = first;
        if (it == last) {
            throw std::invalid_argument{"fromChars called at compile time with empty input"};
        }
        for (; it != last; ++it) {
            if (*it >= '0' && *it <= '9') {
                digit = *it - '0';
            } else if (*it >= 'a' && *it <= 'z') {
                digit = T(10) + (*it - 'a');
            } else if (*it >= 'A' && *it <= 'Z') {
                digit = T(10) + (*it - 'A');
            } else {
                throw std::invalid_argument{"invalid character found"};
            }
            value *= base;
            value += digit;
        }
        return std::from_chars_result{it, {}};
    }
    else {
        return std::from_chars(first, last, value, static_cast<int>(base));
    }
}


constexpr Index roundUpDiv(Index divisor, Index quotient) noexcept {
    assert(divisor >= 0 && quotient > 0);
    return (divisor + quotient - 1) / quotient; // hopefully, this function gets inlined and optimized (quotient is usually a power of 2)
}

ADS_CONSTEVAL static Index bytesNeededForIndexing(Index numElements) noexcept {
    // no constexpr std::bit_floor in C++17; using <= instead of < is fine because no entry actually stores this number
    return numElements <= 256 ? 1 : 1 + bytesNeededForIndexing(roundUpDiv(numElements, 256));
}

template<typename Integer>
ADS_CPP20_CONSTEXPR Integer abs(Integer val) noexcept {
    ADS_IF_CONSTEVAL {
        return val < 0 ? -val : val;
    }
    else {
        return std::abs(val);
    }
}


#ifdef ADS_HAS_CPP20
template<typename Iter1, typename Iter2, typename Cmp = std::compare_three_way>
constexpr auto lexicographicalCompareThreeWay(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, Cmp cmp = Cmp())
        -> decltype(cmp(*first1, *first2)) {
    // libc++ doesn't yet implement lexicographical_compare_three_way
#if !defined _LIBCPP_VERSION
    if constexpr (requires { std::lexicographical_compare_three_way(first1, last1, first2, last2, cmp); }) {
        return std::lexicographical_compare_three_way(first1, last1, first2, last2, cmp);
    }
#else
    for (; first1 != last1 && first2 != last2; ++first1, ++first2) {
        if (auto r = cmp(*first1, *first2); r != 0) {
            return r;
        }
    }
    if (first1 != last1) {
        return std::strong_ordering::greater;
    } else if (first2 != last2) {
        return std::strong_ordering::less;
    } else {
        return std::strong_ordering::equal;
    }
#endif // _LIBCPP_VERSION
}
#endif // ADS_HAS_CPP20


static constexpr Index CACHELINE_SIZE_BYTES = 64;
static constexpr Index ELEMS_PER_CACHELINE = CACHELINE_SIZE_BYTES / 8;


// Useful for testing
inline std::mt19937_64 createRandomEngine() noexcept {
    std::random_device rd;
    std::seed_seq seq{rd(), rd(), rd(), rd()};
    return std::mt19937_64(seq);
}


// Simpler version of std::span, which doesn't exist is C++17
template<typename T>
class Span {
    T* first = nullptr;
    T* last = nullptr;

public:
    using value_type = T;

    constexpr Span() noexcept = default;
    constexpr Span(T* ptr, Index size) noexcept : first(ptr), last(ptr + size) { assert(size >= 0); }
    constexpr Span(T* first, T* last) noexcept : first(first), last(last) { assert(size() >= 0); }

    template<typename Container, typename = std::enable_if_t<std::is_same_v<typename Container::value_type, std::remove_cv_t<T>>>>
    /*implicit*/ ADS_CPP20_CONSTEXPR Span(Container& c) : Span(c.data(), c.size()) {}

    [[nodiscard]] constexpr Index size() const noexcept {
        assert(!std::less<>()(last, first));
        return last - first;
    }

    [[nodiscard]] constexpr T& operator[](Index i) noexcept {
        assert(i >= 0 && i < size());
        return first[i];
    }
    [[nodiscard]] constexpr const T& operator[](Index i) const noexcept {
        assert(i >= 0 && i < size());
        return first[i];
    }
    [[nodiscard]] constexpr bool empty() const noexcept { return size() == 0; }

    constexpr T* data() noexcept { return first; }
    constexpr const T* data() const noexcept { return first; }
    constexpr T* begin() noexcept { return first; }
    constexpr const T* begin() const noexcept { return first; }
    constexpr T* end() noexcept { return last; }
    constexpr const T* end() const noexcept { return last; }

    friend std::ostream& operator<<(std::ostream& os, Span s) noexcept {
        os << "[";
        if (!s.empty()) {
            os << s[0];
        }
        if (s.size() > 0) {
            for (const auto& val : Span(s.begin() + 1, s.end())) {
                os << ", " << val;
            }
        }
        return os;
    }
};

template<typename Container>
Span(const Container&) -> Span<const typename Container::value_type>;

// Simpler version of std::ranges::subrange, C++17 compatible
template<typename Iter, typename Sentinel = Iter>
struct Subrange {
    Iter first;
    Sentinel last;

    using value_type = typename std::iterator_traits<Iter>::value_type;

    [[nodiscard]] constexpr Iter begin() const noexcept { return first; }
    [[nodiscard]] constexpr Sentinel end() const noexcept { return last; }

    [[nodiscard]] constexpr Index size() const noexcept { return last - first; }
};


// Only available in C++20 // TODO: Remove
template<typename T>
std::unique_ptr<T[]> makeUniqueForOverwrite(Index size) noexcept {
    assert(size >= 0);
    if (size == 0) {
        return nullptr;
    }
    return std::unique_ptr<T[]>(new std::remove_extent_t<T>[size]);
}

// TODO: Remove
template<typename T>
ADS_CPP20_CONSTEXPR std::unique_ptr<T[]> toUniquePtr(Span<const T> values, Index size) noexcept {
    std::unique_ptr<T[]> res = makeUniqueForOverwrite<T>(size);
    std::copy(values.begin(), values.end(), res.get());
    // reading uninitialized data would cause UB
    std::fill(res.get() + values.size(), res.get() + size, T());
    return res;
}

template<typename T>
ADS_CPP20_CONSTEXPR std::unique_ptr<T[]> toUniquePtr(Span<const T> values) noexcept {
    return toUniquePtr(values, values.size());
}

} // namespace ads

#endif // ADS_COMMON_HPP
