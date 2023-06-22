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

#ifdef _MSC_VER
#include <intrin.h>
#endif
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
// TODO: Instead of a global Elem alias, define per template to allow smaller sizes
using Elem = std::uint64_t;



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


template<typename T>
ADS_CPP20_CONSTEXPR void uninitializedValueConstructN(T* ADS_RESTRICT dest, Index n) noexcept {
    assert(n >= 0);
    ADS_IF_CONSTEVAL {
        for (Index i = 0; i < n; ++i) {
            std::construct_at(dest + i);
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
constexpr std::from_chars_result fromChars(const char* first, const char* last, T& value, int base = 10) noexcept {
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

template<typename UnsignedInteger> // no concepts in C++17 :(
ADS_CPP20_CONSTEXPR Index log2(UnsignedInteger n) noexcept {
    static_assert(std::is_unsigned_v<UnsignedInteger>);
#ifdef ADS_HAS_CPP20
    return 8 * sizeof(UnsignedInteger) - std::countl_zero(n) - 1;
#elif defined __clang__ || defined __GNUC__
    return 8 * sizeof(unsigned long long) - __builtin_clzll(n) - 1;
#elif defined _MSC_VER
    std::uint64 pos;
    _BitScanReverse(&pos, std::uint64_t(n));
    return 8 * 64 - pos - 1;
#else
    return Index(std::log2(n));
#endif
}

template<typename UnsignedInteger>
ADS_CPP20_CONSTEXPR Index roundUpLog2(UnsignedInteger n) noexcept {
    static_assert(std::is_unsigned_v<UnsignedInteger>);
    assert(n > 0);
#ifdef ADS_HAS_CPP20
    if (std::has_single_bit(n)) {
        return log2(n);
    }
    return log2(n) + 1; // TODO: Test if this is actually faster than the fallback, use for non-c++20 mode as well if faster
#else
    if (n <= 1) {
        return 0;
    }
    return log2(UnsignedInteger(n - 1)) + 1;
#endif
}

template<typename UnsignedInteger>
ADS_CPP20_CONSTEXPR Index popcount(UnsignedInteger n) noexcept {
    static_assert(std::is_unsigned_v<UnsignedInteger>);
#ifdef ADS_HAS_CPP20
    return std::popcount(n);
#elif defined __clang__ || defined __GNUC__
    return __builtin_popcountll(std::uint64_t(n));
#elif defined _MSC_VER
    return __popcnt64(std::uint64_t(n));
#else
    using T = UnsignedInteger;
    // see https://stackoverflow.com/questions/3849337/msvc-equivalent-to-builtin-popcount/42913358#42913358
    n = n - ((n >> 1) & (T) ~(T)0 / 3);                              // temp
    n = (n & (T) ~(T)0 / 15 * 3) + ((n >> 2) & (T) ~(T)0 / 15 * 3);  // temp
    n = (n + (n >> 4)) & (T) ~(T)0 / 255 * 15;                       // temp
    return (T)(n * ((T) ~(T)0 / 255)) >> (sizeof(T) - 1) * CHAR_BIT; // count
#endif
}


template<typename UnsignedInteger>
constexpr UnsignedInteger reverseBits(UnsignedInteger n) noexcept {
    // see https://stackoverflow.com/questions/746171/efficient-algorithm-for-bit-reversal-from-msb-lsb-to-lsb-msb-in-c
    n = ((n & 0xaaaa'aaaa'aaaa'aaaaull) >> 1) | ((n & 0x5555'5555'5555'5555ull) << 1);
    n = ((n & 0xcccc'cccc'cccc'ccccull) >> 2) | ((n & 0x3333'3333'3333'3333ull) << 2);
    if constexpr (sizeof(UnsignedInteger) > 1) {
        n = ((n & 0xf0f0'f0f0'f0f0'f0f0ull) >> 4) | ((n & 0x0f0f'0f0f'0f0f'0f0full) << 4);
    }
    if constexpr (sizeof(UnsignedInteger) > 2) {
        n = ((n & 0xff00'ff00'ff00'ff00ull) >> 8) | ((n & 0x00ff'00ff'00ff'00ffull) << 8);
    }
    if constexpr (sizeof(UnsignedInteger) > 4) {
        n = ((n & 0xffff'0000'ffff'0000ull) >> 16) | ((n & 0x0000'ffff'0000'ffffull) << 16);
    }
    constexpr Index bits = sizeof(UnsignedInteger) * 4;
    return (n >> bits) | (n << bits);
}


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
    /*implicit*/ Span(Container& c) : Span(c.data(), c.size()) { // NOLINT(google-explicit-constructor)
    }

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
    [[nodiscard]] bool empty() const noexcept { return size() == 0; }

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

// Simpler version of std::ranges::subrange
template<typename Iter>
struct Subrange {
    Iter first;
    Iter last;

    using value_type = typename std::iterator_traits<Iter>::value_type;

    [[nodiscard]] Iter begin() const noexcept { return first; }
    [[nodiscard]] Iter end() const noexcept { return last; }

    [[nodiscard]] Index size() const noexcept { return last - first; }
};


// Only available in C++20
template<typename T>
std::unique_ptr<T[]> makeUniqueForOverwrite(Index size) noexcept {
    assert(size >= 0);
    if (size == 0) {
        return nullptr;
    }
    return std::unique_ptr<T[]>(new std::remove_extent_t<T>[size]);
}

template<typename T>
std::unique_ptr<T[]> toUniquePtr(Span<const T> values, Index size) noexcept {
    std::unique_ptr<T[]> res = makeUniqueForOverwrite<T>(size);
    std::copy(values.begin(), values.end(), res.get());
    // reading uninitialized data would cause UB
    std::fill(res.get() + values.size(), res.get() + size, T());
    return res;
}

template<typename T>
std::unique_ptr<T[]> toUniquePtr(Span<const T> values) noexcept {
    return toUniquePtr(values, values.size());
}

} // namespace ads

#endif // ADS_COMMON_HPP
