#ifndef ADS_COMMON_HPP
#define ADS_COMMON_HPP

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
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
#define ADS_ASSUME(x) ADS_ASSUME_IMPL(x)
#else // don't use assume as that could (and probably would) cause the compiler to treat the assert as assert(true)
#define ADS_ASSUME(x) assert(x)
#endif

#ifdef __cpp_lib_assume_aligned // C++20 feature
#define ADS_ASSUME_ALIGNED_IMPL(ptr, align) (ptr) = std::assume_aligned<(align)>(ptr);
#else
#define ADS_ASSUME_ALIGNED_IMPL(ptr, align) /*nothing*/
#endif

#if defined __GNUC__ || defined __clang__
#define ADS_ASSUME_ALIGNED(ptr, align)                                                                                 \
    do {                                                                                                               \
        ADS_IF_CONSTEVAL {}                                                                                            \
        else {                                                                                                         \
            ADS_ASSUME((std::uintptr_t(ptr) & ((align)-1)) == 0);                                                      \
        }                                                                                                              \
        ADS_ASSUME_ALIGNED_IMPL(ptr, align);                                                                           \
    } while (false)
#else
#define ADS_ASSUME_ALIGNED(ptr, align)                                                                                 \
    do {                                                                                                               \
        ADS_IF_CONSTEVAL {}                                                                                            \
        else {                                                                                                         \
            ADS_ASSUME((std::uintptr_t(ptr) & ((align)-1)) == 0);                                                      \
        }                                                                                                              \
        ADS_ASSUME_ALIGNED_IMPL(ptr, align);                                                                           \
    } while (false)
#endif

#ifdef _MSC_VER
#define ADS_FORCE_INLINE(func) __forceinline func
#elif defined __GNUC__ || defined __clang__
#define ADS_FORCE_INLINE(func) inline func __attribute__((always_inline, artificial))
#endif


using Byte = unsigned char;
using Index = std::ptrdiff_t;
using U64 = std::uint64_t;
using I64 = std::int64_t;
using U32 = std::uint32_t;
// TODO: Remove the Elem alias
using Elem = U64;
using Limb = U64;

static constexpr Index CACHELINE_SIZE_BYTES = 64;
static constexpr Index U64_PER_CACHELINE = CACHELINE_SIZE_BYTES / 8;

// Let the compiler assume that accesses are always aligned
template<Index NumBytes = 32>
struct alignas(NumBytes) SIMDLimb {
    constexpr static Index numLimbs = NumBytes / sizeof(Limb);
    static_assert(NumBytes > 0);
    static_assert(NumBytes % sizeof(Limb) == 0);
    static_assert((NumBytes & (NumBytes - 1)) == 0);
    Limb limbs[numLimbs];

    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb& operator[](Index i) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < numLimbs);
        return limbs[i];
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR const Limb& operator[](Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < numLimbs);
        return limbs[i];
    }
};
using CacheLine = SIMDLimb<CACHELINE_SIZE_BYTES>;


struct CreateWithSizeTag {};

struct UninitializedTag {};

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


// inlining this function increases performance in debug mode and helps debugging
[[nodiscard]] ADS_FORCE_INLINE(constexpr Index roundUpDiv(Index dividend, Index divisor) noexcept);
[[nodiscard]] constexpr Index roundUpDiv(Index dividend, Index divisor) noexcept {
    ADS_ASSUME(dividend >= 0 && divisor > 0);
    return (dividend + divisor - 1) / divisor; // hopefully, this function gets inlined and optimized (divisor is usually a power of 2)
}

// inlining this function increases performance in debug mode and helps debugging
[[nodiscard]] ADS_FORCE_INLINE(constexpr Index roundUpTo(Index value, Index divisor) noexcept);
[[nodiscard]] constexpr Index roundUpTo(Index value, Index divisor) noexcept {
    ADS_ASSUME(value >= 0 && divisor > 0);
    return roundUpDiv(value, divisor) * divisor; // hopefully, this function gets inlined and optimized (divisor is usually a power of 2)
}
//
//[[nodiscard]] constexpr Index roundUpTo(Index value, Index divisor1, Index divisor2) noexcept {
//    ADS_ASSUME(value >= 0 && divisor1 > 0 && divisor2 > 0);
//    return roundUpTo(value, std::max(divisor1, divisor2));
//}

[[nodiscard]] ADS_CONSTEVAL static Index bytesNeededForIndexing(Index numElements) noexcept {
    // no constexpr std::bit_floor in C++17; using <= instead of < is fine because no entry actually stores this number
    return numElements <= 256 ? 1 : 1 + bytesNeededForIndexing(roundUpDiv(numElements, 256));
}

template<typename Integer>
[[nodiscard]] ADS_CPP20_CONSTEXPR Integer abs(Integer val) noexcept {
    ADS_IF_CONSTEVAL {
        return val < 0 ? -val : val;
    }
    else {
        return std::abs(val);
    }
}


#ifdef ADS_HAS_CPP20
template<typename Iter1, typename Iter2, typename Cmp = std::compare_three_way>
[[nodiscard]] constexpr auto lexicographicalCompareThreeWay(
        Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, Cmp cmp = Cmp()) -> decltype(cmp(*first1, *first2)) {
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


// Useful for testing
[[nodiscard]] inline std::mt19937_64 createRandomEngine() noexcept {
    std::random_device rd;
    std::seed_seq seq{rd(), rd(), rd(), rd()};
    return std::mt19937_64(seq);
}

namespace detail {
template<Index Val>
struct ConstIndexImpl : TypeIdentity<std::integral_constant<Index, Val>> {};

template<>
struct ConstIndexImpl<(-1)> : TypeIdentity<Index> {};

} // namespace detail

template<Index Val>
using ConstIndex = typename detail::ConstIndexImpl<Val>::Type;

// Simpler version of std::span, which doesn't exist is C++17
template<typename T, Index Size = -1>
class [[nodiscard]] Span {
    using SizeType = ConstIndex<Size>;
    T* first_ = nullptr;
    SizeType size_ = SizeType();

public:
    using value_type = T;

    constexpr Span() noexcept = default;
    constexpr Span(T* ptr, Index size) noexcept : first_(ptr), size_(size) { assert(size >= 0); }
    explicit constexpr Span(T* ptr) noexcept : first_(ptr), size_(SizeType{}) { assert(size() >= 0); }
    template<Index ArrSize>
    explicit constexpr Span(T (&ptr)[ArrSize]) noexcept : first_(ptr), size_(SizeType{}) {
        assert(size() >= 0);
        static_assert(Size <= ArrSize);
    }
    constexpr Span(T* first, T* last) noexcept : first_(first), size_(last - first) { assert(size() >= 0); }

    template<typename Container, typename = std::enable_if_t<std::is_same_v<typename Container::value_type, std::remove_cv_t<T>>>>
    /*implicit*/ ADS_CPP20_CONSTEXPR Span(Container& c) : Span(c.data(), c.size()) {}

    [[nodiscard]] constexpr Index size() const noexcept { return size_; }

    [[nodiscard]] constexpr T& operator[](Index i) noexcept {
        assert(i >= 0 && i < size_);
        return first_[i];
    }
    [[nodiscard]] constexpr const T& operator[](Index i) const noexcept {
        assert(i >= 0 && i < size_);
        return first_[i];
    }
    [[nodiscard]] constexpr bool empty() const noexcept { return size() == 0; }

    [[nodiscard]] constexpr T* data() noexcept { return first_; }
    [[nodiscard]] constexpr const T* data() const noexcept { return first_; }
    [[nodiscard]] constexpr T* begin() noexcept { return first_; }
    [[nodiscard]] constexpr const T* begin() const noexcept { return first_; }
    [[nodiscard]] constexpr T* end() noexcept { return first_ + size_; }
    [[nodiscard]] constexpr const T* end() const noexcept { return first_ + size_; }

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
struct [[nodiscard]] Subrange {
    Iter first;
    Sentinel last;

    using value_type = typename std::iterator_traits<Iter>::value_type;

    [[nodiscard]] constexpr Iter begin() const noexcept { return first; }
    [[nodiscard]] constexpr Sentinel end() const noexcept { return last; }

    [[nodiscard]] constexpr Index size() const noexcept { return last - first; }

    [[nodiscard]] constexpr auto operator[](Index i) const noexcept -> decltype(begin()[i]) { return begin()[i]; }
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
