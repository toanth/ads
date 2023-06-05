#ifndef BITVECTOR_BIT_ACCESS_HPP
#define BITVECTOR_BIT_ACCESS_HPP

#include "common.hpp"

namespace ads {

constexpr static Index dynSize = Index(-1);

template<typename T = Elem, typename Underlying = T>
struct ValuewiseAccess {

    static_assert(sizeof(T) % sizeof(Underlying) == 0);

public:
    [[nodiscard]] static T& getRef(Underlying* ptr, Index i) noexcept {
        if constexpr (std::is_same_v<Underlying, T>) {
            return ptr[i];
        } else {
            return *reinterpret_cast<T*>(ptr);
        }
    }

    [[nodiscard]] static T getBits(const Underlying* ptr, Index i) noexcept {
        assert(i >= 0);
        if constexpr (sizeof(Underlying) == sizeof(T)) {
            return T(ptr[i]);
        } else if constexpr (sizeof(Underlying) * 2 == sizeof(T)) {
            return T(ptr[i]) + (T(ptr[i + 1]) << 8 * sizeof(Underlying));
        } else {
            T val = 0;
            for (Index j = 0; j < sizeof(T) / sizeof(Underlying); ++j) {
                val += T(ptr[i + j]) << (8 * j);
            }
        }
    }

    static void setBits(T* ADS_RESTRICT ptr, Index i, T value) noexcept {
        if constexpr (sizeof(Underlying) == sizeof(T)) {// TODO: Look at code generation to see if these special cases actually help
            ptr[i] = value;
        } else if constexpr (sizeof(Underlying) * 2 == sizeof(T)) {
            ptr[i] = Underlying(value);
            ptr[i + 1] = Underlying(value >> 8);
        } else {
            for (Index j = 0; j < sizeof(T) / sizeof(Underlying); ++j) {
                ptr[i + j] = Underlying(value >> (8 * j));
            }
        }
    }
};

template<Index NumBits, typename T = Elem>
struct BitwiseAccessImpl {

    constexpr static Index bitsInT = sizeof(T) * 8;
    constexpr static T mask = (T(1) << NumBits) - 1;

    static_assert(NumBits > 0 && bitsInT % NumBits == 0);

    [[nodiscard]] static T getBits(const T* ADS_RESTRICT ptr, Index elemStartIndex, Index inElem) noexcept {
        assert(elemStartIndex >= 0 && inElem >= 0 && inElem * NumBits < bitsInT);
        return (ptr[elemStartIndex] >> (NumBits * inElem)) & mask;
    }

    static void setBits(T* ADS_RESTRICT ptr, Index elemStartIndex, Index inElem, T value) noexcept {
        assert(elemStartIndex >= 0 && inElem >= 0 && inElem * NumBits < bitsInT);
        assert(value < (T(1) << NumBits));
        ptr[elemStartIndex] &= ~(mask << (NumBits * inElem));
        ptr[elemStartIndex] |= value << (NumBits * inElem);
    }

public:
    [[nodiscard]] static T getBits(const T* ptr, Index i) noexcept {
        i *= NumBits;
        return getBits(ptr, i / bitsInT, i % bitsInT);
    }

    static void setBits(T* ADS_RESTRICT ptr, Index i, T value) noexcept {
        i *= NumBits;
        setBits(ptr, i / bitsInT, i % bitsInT, value);
    }
};

template<Index NumBits, typename T = Elem>
struct BitwiseAccess : std::conditional_t<NumBits == sizeof(T) * 8, ValuewiseAccess<T>, BitwiseAccessImpl<NumBits, T>> {};

template<typename T>
struct BitwiseAccess<dynSize, T> {
    constexpr static Index bitsInT = sizeof(T) * 8;
    Index numBits = 0;

    [[nodiscard]] Elem getBits(const T* ADS_RESTRICT ptr, Index elemStartIndex, Index bitStartIndex) const noexcept {
        assert(elemStartIndex >= 0 && bitStartIndex >= 0 && bitStartIndex < bitsInT && numBits > 0 && numBits < bitsInT);
        T mask = (T(1) << numBits) - 1;// TODO: Make member?
        T r = (ptr[elemStartIndex] >> bitStartIndex) & mask;
        if (numBits + bitStartIndex > bitsInT) {
            Index remaining = numBits + bitStartIndex - bitsInT;
            mask = (T(1) << remaining) - 1;
            r += (ptr[elemStartIndex + 1] & mask) << (numBits - remaining);
        }
        return r;
    }

    static void setBitsImpl(T* ADS_RESTRICT ptr, Index elemStartIndex, Index bitStartIndex, T value, Index numBits) noexcept {
        assert(numBits > 0 && numBits < bitsInT && bitStartIndex >= 0 && bitStartIndex < bitsInT);
        const T mask = (T(1) << numBits) - 1;
        ptr[elemStartIndex] &= ~(mask << bitStartIndex);
        ptr[elemStartIndex] |= (value & mask) << bitStartIndex;
    }

    void setBits(T* ADS_RESTRICT ptr, Index elemStartIndex, Index bitStartIndex, T value) const noexcept {
        assert(elemStartIndex >= 0 && bitStartIndex >= 0 && numBits > 0 && numBits < bitsInT);
        setBitsImpl(ptr, elemStartIndex, bitStartIndex, value, numBits);
        if (numBits + bitStartIndex > bitsInT) {
            Index remaining = numBits + bitStartIndex - bitsInT;
            setBitsImpl(ptr, elemStartIndex + 1, 0, value >> (numBits - remaining), remaining);
        }
    }

public:
    [[nodiscard]] T getBits(const T* ADS_RESTRICT ptr, Index i) const noexcept {
        i *= numBits;
        return getBits(ptr, i / bitsInT, i % bitsInT);
    }

    void setBits(T* ADS_RESTRICT ptr, Index i, T value) const noexcept {
        i *= numBits;// TODO: Does this get optimized when called in a loop? (ie i += numBits each iteration)
        setBits(ptr, i / bitsInT, i % bitsInT, value);
    }
};

template<Index NumBits = 1, typename T = Elem>
struct BitStorage {
    using BitAccess = BitwiseAccess<NumBits, T>;
    std::unique_ptr<T[]> ptr = nullptr;
    [[no_unique_address]] BitAccess bitAccess = BitAccess();

    [[nodiscard]] Elem getBits(Index i) const noexcept {
        return bitAccess.getBits(ptr.get(), i);
    }

    void setBits(Index i, Elem value) noexcept {
        return bitAccess.setBits(ptr.get(), i, value);
    }

    T operator[](Index i) const noexcept {
        return getBits(i);
    }
};

template<typename T>
using Array = BitStorage<sizeof(T) * 8, T>;


template<Index NumBits = 1, typename T = Elem>
struct BitView {
    using BitAccess = BitwiseAccess<NumBits, T>;
    T* ADS_RESTRICT ptr = nullptr;
    [[no_unique_address]] BitAccess bitAccess = BitAccess();

    BitView() = default;

    template<typename Underlying>
    BitView(Underlying* ptr, Index numT) : ptr(reinterpret_cast<T*>(ptr)) {
        static_assert(alignof(Underlying) % alignof(T) == 0);
        assert(numT >= 0);
        // these two operations should compile to zero assembly instructions for integer types but prevent UB caused by
        // breaking the pointer aliasing rule or by reading/writing to objects whose lifetime hasn't been started yet.
        std::destroy_n(ptr, roundUpDiv(numT, sizeof(Underlying) / sizeof(T)));
        std::uninitialized_value_construct_n(this->ptr, numT);
    }

    [[nodiscard]] T getBits(Index i) const noexcept {
        return bitAccess.getBits(ptr, i);
    }

    void setBits(Index i, T value) noexcept {
        return bitAccess.setBits(ptr, i, value);
    }

    T operator[](Index i) const noexcept {
        return getBits(i);
    }
};

template<typename T>
using View = BitView<sizeof(T) * 8, T>;

}// namespace ads

#endif//BITVECTOR_BIT_ACCESS_HPP
