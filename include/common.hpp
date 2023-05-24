#ifndef ADS_COMMON_HPP
#define ADS_COMMON_HPP

#include <cassert>
#include <cstdint>
#include <cstring>
#include <memory>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace ads {

#if __cplusplus >= 202002L
constexpr static bool hasCpp20 = true;
#define ADS_HAS_CPP20
#define ADS_CONSTEVAL consteval
#define ADS_CPP20_CONSTEXPR constexpr
#include <bit>
#else
constexpr static bool hasCpp20 = false;
#define ADS_CONSTEVAL constexpr
#define ADS_CPP20_CONSTEXPR inline
#endif

using Index = std::ptrdiff_t;

namespace detail {

    template<typename T>
    struct TypeIdentity {// std::type_identity is a C++20 feature
        using Type = T;
    };

    template<Index NumBytes>
    struct IntTypeImpl : TypeIdentity<std::uint64_t> {};

    template<>
    struct IntTypeImpl<4> : TypeIdentity<std::uint32_t> {};
    template<>
    struct IntTypeImpl<3> : TypeIdentity<std::uint32_t> {};
    template<>
    struct IntTypeImpl<2> : TypeIdentity<std::uint16_t> {};
    template<>
    struct IntTypeImpl<1> : TypeIdentity<std::uint8_t> {};

}// namespace detail

template<Index NumBytes>
using IntType = typename detail::IntTypeImpl<NumBytes>::Type;

template<typename Dest>
ADS_CPP20_CONSTEXPR Dest ptrBitCast(const unsigned char* src) noexcept {
#if ADS_HAS_CPP20
    using T = char[sizeof(Dest)];
    return std::bit_cast<Dest>((T*) src);
#else
    Dest res;
    std::memcpy(&res, src, sizeof(Dest));
    return res;
#endif
}

#define ADS_RESTRICT __restrict

template<typename UnsignedInteger>// no concepts in C++17 :(
Index log2(UnsignedInteger n) noexcept {
#ifdef ADS_HAS_CPP20
    return std::bit_floor(n);
#elif defined __clang__ || defined __GNUC__
    return 8 * sizeof(UnsignedInteger) - __builtin_clzll(n) - 1;
#elifdef _MSC_VER
    std::uint64 pos;
    _BitScanReverse(&pos, std::uint64_t(n));
    return 8 * 64 - pos - 1;
#else
    return Index(std::log2(n));
#endif
}

template<typename UnsignedInteger>
Index popcount(UnsignedInteger n) noexcept {
#ifdef ADS_HAS_CPP20
    return std::popcount(n);
#elif defined __clang__ || defined __GNUC__
    return __builtin_popcountll(std::uint64_t(n));
#elifdef _MSC_VER
    return __popcnt64(std::uint64_t(n));
#else
    using T = UnsignedInteger;
    // see https://stackoverflow.com/questions/3849337/msvc-equivalent-to-builtin-popcount/42913358#42913358
    n = n - ((n >> 1) & (T) ~(T) 0 / 3);                              // temp
    n = (n & (T) ~(T) 0 / 15 * 3) + ((n >> 2) & (T) ~(T) 0 / 15 * 3); // temp
    n = (n + (n >> 4)) & (T) ~(T) 0 / 255 * 15;                       // temp
    return (T) (n * ((T) ~(T) 0 / 255)) >> (sizeof(T) - 1) * CHAR_BIT;// count
#endif
}

static constexpr Index CACHELINE_SIZE_BYTES = 64;
static constexpr Index ELEMS_PER_CACHELINE = CACHELINE_SIZE_BYTES / 8;

// The following functions are only available in C++20 :(
template<typename T>
std::unique_ptr<T[]> makeUniqueForOverwrite(Index size) noexcept {
    return std::unique_ptr<T[]>(new std::remove_extent_t<T>[size]);
}

constexpr Index roundUpDiv(Index divisor, Index quotient) noexcept {
    return (divisor + quotient - 1) / quotient;// hopefully, this function gets inlined and optimized (quotient is usually a power of 2)
}

// Simplified version of std::span
template<typename T>
class Span {
    T* first = nullptr;
    T* last = nullptr;

public:
    constexpr Span() noexcept = default;
    constexpr Span(T* ptr, Index size) noexcept : first(ptr), last(ptr + size) {
        assert(size >= 0);
    }
    constexpr Span(T* first, T* last) noexcept : first(first), last(last) {
        assert(size() >= 0);
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

    constexpr T* begin() noexcept { return first; }
    constexpr const T* begin() const noexcept { return first; }
    constexpr T* end() noexcept { return last; }
    constexpr const T* end() const noexcept { return last; }
};

}// namespace ads

#endif//ADS_COMMON_HPP
