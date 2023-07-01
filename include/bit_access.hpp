#ifndef BITVECTOR_BIT_ACCESS_HPP
#define BITVECTOR_BIT_ACCESS_HPP

#include "common.hpp"
#include <utility>

namespace ads {

constexpr static Index dynSize = Index(-1);

//---- Templates that represent various ways to access sequences of bits in an array ----

/// \brief Access entire integer-like objects of type T, from an array of type Underlying[]
/// \tparam T The logical type of array elements
/// \tparam Underlying The actual type of array elements as they are stored (which may again be different from the type that was allocated)
template<typename T = Elem, typename Underlying = T>
struct ValuewiseAccess {

    static_assert(sizeof(T) % sizeof(Underlying) == 0);

public:
    [[nodiscard]] constexpr static T& getRef(Underlying* ptr, Index i) noexcept {
        if constexpr (std::is_same_v<Underlying, T>) {
            return ptr[i];
        } else {
            return *reinterpret_cast<T*>(ptr);
        }
    }

    [[nodiscard]] static constexpr T getBits(const Underlying* ptr, Index i) noexcept {
        assert(i >= 0);
        assert(ptr);
        if constexpr (sizeof(Underlying) == sizeof(T)) {
            return T(ptr[i]);
        } else {
            T val = 0;
            for (Index j = 0; j < sizeof(T) / sizeof(Underlying); ++j) {
                val += T(ptr[i + j]) << (j * 8 * sizeof(Underlying));
            }
            return val;
        }
    }

    static constexpr void setBits(T* ADS_RESTRICT ptr, Index i, T value) noexcept {
        assert(i >= 0);
        assert(ptr);
        if constexpr (sizeof(Underlying) == sizeof(T)) { // TODO: Look at code generation to see if these special cases actually help
            ptr[i] = value;
        } else {
            for (Index j = 0; j < sizeof(T) / sizeof(Underlying); ++j) {
                ptr[i + j] = Underlying(value >> (j * 8 * sizeof(Underlying)));
            }
        }
    }
};


/// \brief Access N-bit sequences
/// \tparam NumBits The number of bits in a single element, not necessarily a multiple of 8
/// \tparam T The type of the underling array elements
template<Index NumBits, typename T = Elem>
struct BitwiseAccessImpl {

    constexpr static Index bitsInT = sizeof(T) * 8;
    constexpr static T mask = (T(1) << NumBits) - 1;

    static_assert(NumBits > 0 && bitsInT % NumBits == 0);

    [[nodiscard]] static constexpr T getBits(const T* ADS_RESTRICT ptr, Index elemStartIndex, Index inElem) noexcept {
        assert(elemStartIndex >= 0 && inElem >= 0 && inElem * NumBits < bitsInT);
        assert(ptr);
        return (ptr[elemStartIndex] >> (NumBits * inElem)) & mask;
    }

    static constexpr void setBits(T* ADS_RESTRICT ptr, Index elemStartIndex, Index inElem, T value) noexcept {
        assert(elemStartIndex >= 0 && inElem >= 0 && inElem * NumBits < bitsInT);
        assert(value < (T(1) << NumBits));
        assert(ptr);
        ptr[elemStartIndex] &= ~(mask << (NumBits * inElem));
        ptr[elemStartIndex] |= value << (NumBits * inElem);
    }

public:
    [[nodiscard]] static constexpr T getBits(const T* ptr, Index i) noexcept {
        i *= NumBits;
        return getBits(ptr, i / bitsInT, i % bitsInT);
    }

    static constexpr void setBits(T* ADS_RESTRICT ptr, Index i, T value) noexcept {
        i *= NumBits;
        setBits(ptr, i / bitsInT, i % bitsInT, value);
    }
};


template<Index NumBits, typename T = Elem>
struct BitwiseAccess
    : std::conditional_t<NumBits % (sizeof(T) * 8) == 0, ValuewiseAccess<T>, BitwiseAccessImpl<NumBits, T>> {};


/// \brief Access N-bit sequences where n is only known at runtime
/// \tparam T The type of the underlying array elements.
template<typename T>
struct BitwiseAccess<dynSize, T> {
    constexpr static Index bitsInT = sizeof(T) * 8;
    Index numBits = 0;

    [[nodiscard]] constexpr Elem getBits(const T* ADS_RESTRICT ptr, Index elemStartIndex, Index bitStartIndex) const noexcept {
        assert(elemStartIndex >= 0 && bitStartIndex >= 0 && bitStartIndex < bitsInT && numBits > 0 && numBits < bitsInT);
        assert(ptr);
        T mask = (T(1) << numBits) - 1; // TODO: Make member?
        T r = (ptr[elemStartIndex] >> bitStartIndex) & mask;
        if (numBits + bitStartIndex > bitsInT) {
            Index remaining = numBits + bitStartIndex - bitsInT;
            mask = (T(1) << remaining) - 1;
            r += (ptr[elemStartIndex + 1] & mask) << (numBits - remaining);
        }
        return r;
    }

    static constexpr void setBitsImpl(T* ADS_RESTRICT ptr, Index elemStartIndex, Index bitStartIndex, T value, Index numBits) noexcept {
        assert(numBits > 0 && numBits < bitsInT && bitStartIndex >= 0 && bitStartIndex < bitsInT);
        const T mask = (T(1) << numBits) - 1;
        ptr[elemStartIndex] &= ~(mask << bitStartIndex);
        ptr[elemStartIndex] |= (value & mask) << bitStartIndex;
    }

    void constexpr setBits(T* ADS_RESTRICT ptr, Index elemStartIndex, Index bitStartIndex, T value) const noexcept {
        assert(elemStartIndex >= 0 && bitStartIndex >= 0 && numBits > 0 && numBits < bitsInT);
        assert(ptr);
        setBitsImpl(ptr, elemStartIndex, bitStartIndex, value, numBits);
        if (numBits + bitStartIndex > bitsInT) {
            Index remaining = numBits + bitStartIndex - bitsInT;
            setBitsImpl(ptr, elemStartIndex + 1, 0, value >> (numBits - remaining), remaining);
        }
    }

public:
    [[nodiscard]] constexpr T getBits(const T* ADS_RESTRICT ptr, Index i) const noexcept {
        i *= numBits;
        return getBits(ptr, i / bitsInT, i % bitsInT);
    }

    constexpr void setBits(T* ADS_RESTRICT ptr, Index i, T value) const noexcept {
        i *= numBits; // TODO: Does this get optimized when called in a loop? (ie i += numBitsPerValue each iteration)
        setBits(ptr, i / bitsInT, i % bitsInT, value);
    }
};

/// \brief Allocate uninitialized memory
/// \tparam T element type
/// \param n number of elements
/// \return uninitialized memory, represented as an array of T with n elements
// although allocation isn't technically noexcept, there's no reason not to crash the program in that case
template<typename T>
ADS_CPP20_CONSTEXPR T* allocate(Index n) noexcept {
    std::allocator<T> alloc;
    return alloc.allocate(n); // contexpr since C++20
}

/// \brief Deallocate memory previously allocated by allocate()
/// \tparam T element type
/// \param ptr allocated array
/// \param n number of elements
template<typename T>
ADS_CPP20_CONSTEXPR void deallocate(T* ptr, Index n) noexcept {
    std::allocator<T> alloc;
    if (ptr) {
        alloc.deallocate(ptr, n); // contexpr since C++20
    }
}


//---- Tying it all together: The BitView template and its variations combine storing an array with accessing bits ----


/// \brief Class template which represents an array of N-bit values and manages its elements' lifetimes, but doesn't
/// concern itself with the (de)allocation of its memory; can be constructed from an uninitialized array of another
/// integer type.
/// \tparam NumBits The number of bits per element, doesn't have to be a multiple of 8. -1 is used if the size is only known at runtime.
/// \tparam T The type of a single array element. Often, NumBits is sizeof(T) * 8, in which case accessing an element is simply an array access.
template<Index NumBits = 1, typename T = Elem>
struct BitView {
    using BitAccess = BitwiseAccess<NumBits, T>;
    T* ADS_RESTRICT ptr = nullptr; // Not a unique_ptr because this may only be one part of an allocation
    Index numT = 0;
    [[no_unique_address]] BitAccess bitAccess = BitAccess();

    constexpr BitView() noexcept = default;

    ADS_CPP20_CONSTEXPR ~BitView() noexcept {
        std::destroy_n(ptr, numT); // should compile to 0 instructions for integer types
        ADS_IF_CONSTEVAL {
            deallocate(ptr, numT);
            ptr = nullptr;
        }
    }

    /// \brief Construct the BitView from a pointer to *uninitialized* memory (such as from allocate()) of type Underlying.
    /// \tparam Underlying A pointer to Underlying must be convertible to a pointer to T.
    /// \param underlyingPtr will be converted to a T* and used as the begin of the BitView.
    /// \param numT the number of array elements of type T in the BitView. Not necessarily the same as the original array's size.
    template<typename Underlying>
    ADS_CPP20_CONSTEXPR BitView(Underlying* underlyingPtr, Index numT) : numT(numT) {
        static_assert(alignof(Underlying) % alignof(T) == 0);
        ADS_IF_CONSTEVAL { // reinterpret_cast isn't constexpr until probably C++26
            // always reallocate because else it would get tricky to figure out whether to deallocate in the destructor
            ptr = allocate<T>(numT); // memory usage or performance isn't a big concern for constant evaluation
        }
        else if constexpr (std::is_same_v<Underlying, T>) {
            ptr = underlyingPtr;
        }
        else {
            ptr = reinterpret_cast<T*>(underlyingPtr);
        }
        assert(numT >= 0);
        // the following call should compile to zero assembly instructions for integer types but prevent UB caused by
        // breaking the pointer aliasing rule or by reading/writing to objects whose lifetime hasn't been started yet.
        uninitializedValueConstructN(ptr, numT);
    }

    ADS_CPP20_CONSTEXPR BitView(BitView&& other) noexcept
        : ptr(std::exchange(other.ptr, nullptr)), numT(std::exchange(other.numT, 0)), bitAccess(other.bitAccess) {}

    ADS_CPP20_CONSTEXPR BitView& operator=(BitView&& other) noexcept {
        using std::swap;
        swap(ptr, other.ptr);
        swap(numT, other.numT);
        bitAccess = other.bitAccess;
        return *this;
    }

    [[nodiscard]] constexpr Index numBitsPerValue() const noexcept {
        if constexpr (NumBits == -1) {
            return bitAccess.numBits;
        } else {
            return NumBits;
        }
    }

    [[nodiscard]] constexpr T getBits(Index i) const noexcept {
        assert(i >= 0);
        assert((i + 1) * numBitsPerValue() <= numT * sizeof(T) * 8);
        return bitAccess.getBits(ptr, i);
    }

    constexpr void setBits(Index i, T value) noexcept {
        assert(i >= 0);
        assert((i + 1) * numBitsPerValue() <= numT * sizeof(T) * 8);
        return bitAccess.setBits(ptr, i, value);
    }

    constexpr T operator[](Index i) const noexcept { return getBits(i); }
};

template<typename T>
using View = BitView<sizeof(T) * 8, T>;

template<typename T = Elem>
class Allocation {
    T* mem = nullptr;
    Index numTs = 0;

public:
    ADS_CPP20_CONSTEXPR Allocation() noexcept = default;
    ADS_CPP20_CONSTEXPR Allocation(Index n, T* ptr = nullptr) noexcept : numTs(n) {
        assert(n >= 0);
        if (!ptr) {
            mem = allocate<T>(n);
        } else {
            mem = ptr;
            numTs = -numTs;
        }
    }

    ADS_CPP20_CONSTEXPR Allocation(Allocation&& other) noexcept
        : mem(std::exchange(other.mem, nullptr)), numTs(std::exchange(other.numTs, 0)) {}

    ADS_CPP20_CONSTEXPR ~Allocation() noexcept {
        if (isOwner()) {
            deallocate(memory(), size());
        }
    }

    ADS_CPP20_CONSTEXPR Allocation& operator=(Allocation&& other) noexcept {
        using std::swap;
        swap(mem, other.mem);
        swap(numTs, other.numTs);
        return *this;
    }

    [[nodiscard]] constexpr Index size() const noexcept { return abs(numTs); }

    [[nodiscard]] constexpr T* memory() const noexcept { return mem; }

    [[nodiscard]] constexpr bool isOwner() const noexcept { return numTs > 0; }
};


///// Like BitView<T>, but manages (de)allocation of its memory, similar to a combination of Allocation<T> and BitView<T>.
// template<typename T, Index NumBits = 8 * sizeof(T)>
// class Array : public BitView<NumBits, T> {
//     using Base = BitView<NumBits, T>;
//
//     constexpr Array(CreateWithSizeTag, Index numTs) : Base(allocate<T>(numTs), numTs) {}
//
// public:
//     constexpr Array() noexcept = default;
//
//     explicit constexpr Array(Index n) noexcept : Array(CreateWithSizeTag(), roundUpDiv(n * NumBits, (8 * sizeof(T)))) {}
//
//     constexpr Array(Index n, Index bitLen) noexcept
//         : Array(CreateWithSizeTag(), roundUpDiv(n * bitLen, (8 * sizeof(T)))) {
//         static_assert(NumBits == dynSize, "Don't use this constructor for fixed-sized bit views");
//         this->bitAccess.numBitsPerValue = bitLen;
//     }
//
//     ADS_CPP20_CONSTEXPR ~Array() noexcept { deallocate(this->ptr, this->numT); }
// };
//
// template<Index NumBits>
// using BitStorage = Array<Elem, NumBits>;

} // namespace ads

#endif // BITVECTOR_BIT_ACCESS_HPP
