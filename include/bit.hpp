#ifndef ADS_BIT_HPP
#define ADS_BIT_HPP

#include "common.hpp"
#include <array>


namespace ads {

template<typename Integer> // no concepts in C++17 :(
[[nodiscard]] ADS_CPP20_CONSTEXPR Index intLog2(Integer n) noexcept {
    if constexpr (std::is_signed_v<Integer>) {
        // converting a negative value to an unsigned integer was implementation defined before C++20, but did the right thing
        return intLog2(static_cast<std::make_unsigned_t<Integer>>(n));
    } else {
        static_assert(std::is_unsigned_v<Integer>);
        ADS_ASSUME(n > 0); // helps to eliminate a check in std::countl_zero; __builtin_clz produces UB in that case
#ifdef ADS_HAS_CPP20
        return 8 * sizeof(Integer) - std::countl_zero(n) - 1;
#elif defined ADS_HAS_DEFAULT_GCC_INTRINSICS
        return 8 * sizeof(unsigned long long) - __builtin_clzll(n) - 1;
#elif defined ADS_HAS_MSVC_INTRINSICS
        std::uint64 pos;
        _BitScanReverse64(&pos, std::uint64_t(n));
        return 64 - pos - 1;
#else
        return Index(std::log2(n));
#endif
    }
}

template<typename Integer>
[[nodiscard]] ADS_CPP20_CONSTEXPR Index roundUpLog2(Integer n) noexcept {
    if constexpr (std::is_signed_v<Integer>) {
        return roundUpLog2(static_cast<std::make_unsigned_t<Integer>>(n));
    } else {
        static_assert(std::is_unsigned_v<Integer>);
        assert(n > 0);
#ifdef ADS_HAS_CPP20
        if (std::has_single_bit(n)) {
            return intLog2(n);
        }
        return intLog2(n) + 1; // TODO: Test if this is actually faster than the fallback, use for non-c++20 mode as well if faster
#else
        if (n <= 1) {
            return 0;
        }
        return intLog2(Integer(n - 1)) + 1;
#endif
    }
}

template<typename T>
[[nodiscard]] constexpr Index popcountFallback(T n) noexcept {
    static_assert(std::is_unsigned_v<T>);
    // see https://stackoverflow.com/questions/3849337/msvc-equivalent-to-builtin-popcount/42913358#42913358,
    // or https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
    n = n - ((n >> 1) & (T) ~(T)0 / 3);                             // temp
    n = (n & (T) ~(T)0 / 15 * 3) + ((n >> 2) & (T) ~(T)0 / 15 * 3); // temp
    n = (n + (n >> 4)) & (T) ~(T)0 / 255 * 15;                      // temp
    return (T)(n * ((T) ~(T)0 / 255)) >> (sizeof(T) - 1) * 8;       // count
}

template<typename UnsignedInteger>
[[nodiscard]] ADS_CPP20_CONSTEXPR Index popcount(UnsignedInteger n) noexcept {
    static_assert(std::is_unsigned_v<UnsignedInteger>);
#ifdef ADS_HAS_CPP20
    return std::popcount(n);
#elif defined ADS_HAS_DEFAULT_GCC_INTRINSICS
    return __builtin_popcountll(std::uint64_t(n));
#elif defined ADS_HAS_MSVC_INTRINSICS
    return __popcnt64(std::uint64_t(n));
#else
    return popcountFallback(n);
#endif
}

/// \brief Returns the rank one of pos in val, ie. counts all ones in (pos, lsb]
template<typename UnsignedInteger>
[[nodiscard]] ADS_CPP20_CONSTEXPR Index popcountBefore(UnsignedInteger val, Index pos) noexcept {
    ADS_ASSUME(pos >= 0);
    ADS_ASSUME(pos < Index(sizeof(UnsignedInteger) * 8));
    const Index shiftAmount = sizeof(UnsignedInteger) * 8 - pos - 1;
    return popcount((val << shiftAmount) << 1); // two shifts to ignore bit pos without invoking UB
}

template<typename UnsignedInteger>
[[nodiscard]] constexpr UnsignedInteger reverseBits(UnsignedInteger n) noexcept {
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


[[nodiscard]] ADS_CPP20_CONSTEXPR Index countTrailingZeros(U64 n) noexcept {
    ADS_ASSUME(n > 0); // undefined otherwise
#ifdef ADS_HAS_CPP20
    return static_cast<Index>(std::countr_zero(n));
#elif defined ADS_HAS_DEFAULT_GCC_INTRINSICS
    return static_cast<Index>(__builtin_ctzll(n));
    return static_cast<Index>(__tzcnt_u64(n));
#elif defined ADS_HAS_MSVC_INTRINSICS
    return static_cast<Index>(_BitScanForward64(n));
#else
    return static_cast<Index>(64 - 1 - intLog2(reverseBits(n))); // not very fast, but this shouldn't really be executed anyway
#endif
}

[[nodiscard]] ADS_CPP20_CONSTEXPR U64 u64SelectImpl(U64 n, Index bitRank) noexcept {
    // #if defined ADS_HAS_GCC_BMI2 || (defined ADS_HAS_MSVC_INTRINSICS && defined ADS_USE_BMI2_INTRINSICS)
    //     // see https://stackoverflow.com/questions/7669057/find-nth-set-bit-in-an-int/27453505#27453505
    //     // TODO: Measure if faster than fallback (in general, but especially) on AMD processors before Zen 3.
    //     // At least on my machine, uncommenting this results in a ridiculously slow (up to 40x slower for small bvs) select
    //     // because the processor advertises BMI2 as supported but implements it in microcode
    //     // It's surprisingly hard to figure out at compile time whether the target architecture supports BMI2 instructions,
    //     // so let the user decide (with the default being the generic fallback) in the MSVC case
    //     return _pdep_u64(1ull << bitRank, n);
    // #else
    //  TODO: Use lookup table? Measure!
    //  see https://stackoverflow.com/questions/7669057/find-nth-set-bit-in-an-int/7669326#7669326
    for (Index i = 0; i < bitRank; ++i) {
        n &= n - 1;
    }
    return n & ~(n - 1);
    // #endif
}


// Depending on macros, not every version of u64Select is constexpr and there is no std::is_constant_evaluated in C++17
[[nodiscard]] ADS_CONSTEVAL Index constevalU64Select(U64 value, Index bitRank) noexcept {
    Index numSet = popcountFallback(value);
    if (bitRank >= numSet) {
        return -1;
    }
    for (Index i = 0; i < 64; ++i) {
        value >>= 1;
        if (numSet - popcountFallback(value) > bitRank) {
            return i;
        }
    }
    assert(false);
    return -1;
}

template<Index BitSize = 8>
using BitSelectTable = std::array<std::array<Byte, BitSize>, (1ull << BitSize)>;

template<Index BitSize = 8>
[[nodiscard]] ADS_CONSTEVAL BitSelectTable<BitSize> precomputeBitSelectTable() noexcept {
    BitSelectTable<BitSize> table{}; // {} needed for constant evaluation
    for (Index bitString = 0; bitString < Index(table.size()); ++bitString) {
        for (Index bitRank = 0; bitRank < BitSize; ++bitRank) {
            table[bitString][bitRank] = static_cast<Byte>(constevalU64Select(bitString, bitRank));
        }
    }
    return table;
}

constexpr static inline BitSelectTable<> byteSelectTable = precomputeBitSelectTable();

[[nodiscard]] constexpr Index byteSelectWithTable(Byte byte, Index bitRank) noexcept {
    return byteSelectTable[byte][bitRank];
}

[[nodiscard]] ADS_CPP20_CONSTEXPR Index u64SelectWithTable(U64 value, Index bitRank) noexcept {
    // code based on https://github.com/s-yata/marisa-trie/blob/master/lib/marisa/grimoire/vector/bit-vector.cc#L180
    // this effectively performs a parallel popcount (compare https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
    // or the popcountFallback() implementation) for each individual byte in parallel,
    // then selects the correct byte and uses the lookup table to select within the chosen byte.
    U64 counts = value - ((value >> 1) & 0x5555'5555'5555'5555ull);
    counts = (counts & 0x3333'3333'3333'3333ull) + ((counts >> 2) & 0x3333'3333'3333'3333ull);
    counts = (counts + (counts >> 4)) & 0x0f0f'0f0f'0f0f'0f0full;
    counts *= 0x0101'0101'0101'0101ull; // parallel prefix sum
    // Now, find the first byte with higher popcount prefix sum: The maximum (rank + 1) is 64, so add 128 to each byte
    // and then subtract the bytewise prefix count. The rightmost byte where the 7th bit remained set is the first byte
    // where bitRank + 1 is greater than the number of set bits up to and including this byte.
    U64 x = (counts | 0x8080'8080'8080'8080ull) - ((bitRank + 1) * 0x0101'0101'0101'0101ull);
    Index numBitsBeforeByte = countTrailingZeros((x & 0x8080'8080'8080'8080ull) >> 7);
    ADS_ASSUME(numBitsBeforeByte >= 0);
    ADS_ASSUME(numBitsBeforeByte <= 64 - 8);
    bitRank -= Index(((counts << 8) >> numBitsBeforeByte) & 0xff);
    ADS_ASSUME(bitRank >= 0);
    ADS_ASSUME(bitRank < 8);
    return numBitsBeforeByte + byteSelectWithTable(static_cast<Byte>(value >> numBitsBeforeByte), bitRank);
}


[[nodiscard]] ADS_CPP20_CONSTEXPR Index u64Select(U64 value, Index bitRank) noexcept {
    return u64SelectWithTable(value, bitRank);
    //    return countTrailingZeros(u64SelectImpl(value, bitRank));
}


[[nodiscard]] ADS_CPP20_CONSTEXPR Index u256Select(const U64* valuePtr, Index bitRank) noexcept {
    ADS_ASSUME(bitRank >= 0);
    ADS_ASSUME(bitRank < 256);
    ADS_ASSUME(valuePtr);
    ADS_ASSUME_ALIGNED(valuePtr, 32);
    Index count = popcount(valuePtr[0]);
    if (bitRank < count) {
        return u64Select(valuePtr[0], bitRank);
    }
    bitRank -= count;
    count = popcount(valuePtr[1]);
    if (bitRank < count) {
        return 64 + u64Select(valuePtr[1], bitRank);
    }
    bitRank -= count;
    count = popcount(valuePtr[2]);
    if (bitRank < count) {
        return 128 + u64Select(valuePtr[2], bitRank);
    }
    bitRank -= count;
    ADS_ASSUME(bitRank >= 0);
    ADS_ASSUME(bitRank < 64);
    ADS_ASSUME(bitRank < popcount(valuePtr[3]));
    return 3 * 64 + u64Select(valuePtr[3], bitRank);
}

[[nodiscard]] ADS_CPP20_CONSTEXPR Index u256SelectZero(const U64* valuePtr, Index bitRank) noexcept {
    ADS_ASSUME(bitRank >= 0);
    ADS_ASSUME(bitRank < 256);
    ADS_ASSUME(valuePtr);
    ADS_ASSUME_ALIGNED(valuePtr, 32);
    alignas(32) U64 negated[] = {~valuePtr[0], ~valuePtr[1], ~valuePtr[2], ~valuePtr[3]};
    return u256Select(negated, bitRank);
}

[[nodiscard]] ADS_CPP20_CONSTEXPR Index u256Rank(const U64* valuePtr) noexcept {
    ADS_ASSUME(valuePtr);
    ADS_ASSUME_ALIGNED(valuePtr, 32);
    // TODO: Use SSE instructions?
    return popcount(valuePtr[0]) + popcount(valuePtr[1]) + popcount(valuePtr[2]) + popcount(valuePtr[3]);
}

[[nodiscard]] ADS_CPP20_CONSTEXPR Index u256RankBefore(const U64* valuePtr, Index pos) noexcept {
    ADS_ASSUME(valuePtr);
    ADS_ASSUME_ALIGNED(valuePtr, 32);
    ADS_ASSUME(pos >= 0);
    ADS_ASSUME(pos < 256);
    Index res = 0;
    for (Index i = 0; i < pos / 64; ++i) {
        res += popcount(valuePtr[i]);
    }
    return res + popcountBefore(valuePtr[pos / 64], pos % 64);
}



} // namespace ads

#endif // ADS_BIT_HPP
