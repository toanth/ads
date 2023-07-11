#ifndef ADS_BITVEC_HPP
#define ADS_BITVEC_HPP

#include "concept.hpp"

namespace ads {


/// \brief CRTP base class of all bitvectors. Most bitvectors inherit from NormalBitvecBase instead, which inherits
/// from this template. The CRTP works like "normal" inheritance with virtual methods, except that "virtual
/// dispatch" 1. is opt-in by calling derived().func(), 2. works for static methods and templates, 3. is resolved at
/// compile time, and 4. allows signatures to differ. On the downside, there is no single baseclass, so there is no
/// built-in type erasure over different bitvector implementations. The main benefit over simply using a concept is
/// that this allows default implementations for all methods. The `derived()` calls are admittedly ugly but would be
/// unnecessary with C++23's deduced this. \tparam Derived The (most) derived Bitvector class, which needs to
/// inherit from this template like this: `class ActualBitvector : public BitvecBase<ActualBitvector>`. Inheritance
/// can be private; the base should be declared as a friend of the derived class.
template<typename Derived>
class [[nodiscard]] BitvecBase {
protected:
    friend Derived;
    /// This owns the uninitialized allocated memory. Being in the base class ensures the memory is (de)allocated at the right time.
    /// In many cases, such as for the Elias-Fano class, this will hold nullptr because the actual owner is another object.
    Allocation<> allocation = Allocation();

    Index numBits_ = 0;

    // ** The derived() methods allow simulation of virtual dispatch at compile time **
    // Inlining these functions even in debug builds improves speed and debugging experience
    [[nodiscard]] ADS_FORCE_INLINE(ADS_CPP20_CONSTEXPR Derived& derived() noexcept) {
        return static_cast<Derived&>(*this);
    }
    [[nodiscard]] ADS_FORCE_INLINE(ADS_CPP20_CONSTEXPR const Derived& derived() const noexcept) {
        return static_cast<const Derived&>(*this);
    }

    // ** Construction from string_view **
    ADS_CPP20_CONSTEXPR static auto limbViewFromStringView(const std::string_view& str, Index base) {
        if (base != 2 && base != 4 && base != 16) [[unlikely]] {
            throw std::invalid_argument("base must be one of 2, 4, or 16");
        }
        const Index log2ofBase = intLog2(base);
        const Index charsPerLimb = 64 / log2ofBase;
        auto proj = [base, log2ofBase, charsPerLimb](std::string_view str, Index i) -> Limb {
            ADS_ASSUME(str.size() > i * charsPerLimb);
            str.remove_prefix(i * charsPerLimb);
            Index numToParse = std::min(charsPerLimb, Index(str.size()));
            U64 res;
            auto err = fromChars(str.data(), str.data() + numToParse, res, int(base));
            if (err.ec != std::errc()) [[unlikely]] {
                throw std::invalid_argument(std::make_error_code(err.ec).message());
            } else if (err.ptr != str.data() + numToParse) [[unlikely]] {
                throw std::invalid_argument("invalid character found");
            }
            if (numToParse < charsPerLimb) {
                res <<= (charsPerLimb - numToParse) * log2ofBase;
            }
            return reverseBits(res);
        };
        using Iter = RandAccessIter<std::string_view, decltype(proj)>;
        auto begin = Iter(str, proj);
        auto end = Iter(str, proj, roundUpDiv(Index(str.size()), charsPerLimb));
        return Subrange<Iter>(begin, end);
    }

    /// \brief Allocates memory for the Bitvector or uses the provided ptr, assuming it points to enough memory.
    /// \param numBits The number of bits stored in this bitvector, usually less than the total allocated space.
    /// \param ptr If nullptr, allocate memory. Else, use this pointer; it must point to enough memory.
    template<typename Underlying = CacheLine>
    ADS_CPP20_CONSTEXPR BitvecBase(Index numBits, Underlying* ptr = nullptr) noexcept
        : allocation(Derived::template allocatedSizeForBitsIn<Group::CacheLine>(numBits) * CACHELINE_SIZE_BYTES, ptr),
          numBits_(numBits) {}

public:
    constexpr BitvecBase() noexcept = default;

    static ADS_CPP20_CONSTEXPR Derived uninitializedForSize(Index numBits, CacheLine* mem = nullptr) noexcept {
        return Derived(UninitializedTag{}, numBits, mem);
    }

    template<Group G>
    [[nodiscard]] static constexpr Index allocatedSizeForBitsIn(Index numBits) noexcept {
        Index inBits = Derived::allocatedSizeInBytesForBits(numBits) * 8;
        if constexpr (G == Group::CacheLine) {
            return roundUpDiv(inBits, CACHELINE_SIZE_BYTES * 8);
        } else {
            return roundUpDiv(inBits, Derived::template numBitsIn<G>());
        }
    }

    [[nodiscard]] const Allocation<>& alloc() const noexcept { return allocation; }


    // must be implemented in the derived class
    //    [[nodiscard]] static constexpr Index allocatedSizeInBytesForBits(Index numBits) noexcept {}

    template<Group G>
    [[nodiscard]] constexpr Index allocatedSizeIn() const noexcept {
        Index inBits = derived().allocatedSizeInBytes() * 8;
        ADS_ASSUME(allocation.sizeInBytes() * 8 == inBits);
        if constexpr (G == Group::CacheLine) {
            return roundUpDiv(inBits, CACHELINE_SIZE_BYTES * 8);
        } else {
            return roundUpDiv(inBits, Derived::template numBitsIn<G>());
        }
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index allocatedSizeInBytes() const noexcept {
        Index inBytes = allocation.sizeInBytes();
        ADS_ASSUME(inBytes == Derived::allocatedSizeInBytesForBits(derived().size()));
        ADS_ASSUME(inBytes % Derived::requiredAlignment() == 0);
        return inBytes;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index allocatedSizeInBits() const noexcept {
        return derived().allocatedSizeInBytes() * 8;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return numBits_; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBits() const noexcept { return numBits_; }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numOnes() const noexcept {
        return derived().numBits() - derived().numZeros();
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numZeros() const noexcept {
        return derived().numBits() - derived().numOnes();
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBitsEqualTo() const noexcept {
        if constexpr (IsOne) {
            return derived().numOnes();
        } else {
            return derived().numZeros();
        }
    }

    // The following three templated functions shouldn't be overwritten by a derived class (except for SuperblockBitvec,
    // which adds a case for Group::Superblock to numBitsIn)
    template<Group G>
    [[nodiscard]] static constexpr Index numBitsIn() noexcept {
        if constexpr (G == Group::Bit) {
            return 1;
        } else if constexpr (G == Group::Byte) {
            return 8;
        } else if constexpr (G == Group::Limb) {
            return 64;
        } else if constexpr (G == Group::Block) {
            return Derived::blockSize();      // compile error for non-NormalBitvecs
        } else if constexpr (G == Group::CacheLine) {
            return Derived::cacheLineSize();  // compile error for non-NormalBitvecs
        } else if constexpr (G == Group::Superblock) {
            return Derived::superblockSize(); // compile error for non-SuperblockBitvecs
        } else {
            static_assert(G != G, "Unknown group size");
        }
    }

    template<Group G>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index sizeIn() const noexcept {
        return roundUpDiv(derived().size(), numBitsIn<G>());
    }

    template<Group G>
    [[nodiscard]] constexpr bool get(Index bitIdx) const noexcept {
        static_assert(G == Group::Bit, "Only a NormalBitvec uses other groups than bits");
        return derived().getBit(bitIdx);
    }

    using BitIter = RandAccessIter<Derived, decltype(getGroupFunc<Group::Bit>)>;

    [[nodiscard]] ADS_CPP20_CONSTEXPR BitIter bitIter(Index i) const {
        return BitIter(derived(), getGroupFunc<Group::Bit>, i);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<BitIter> bitView() const noexcept {
        return Subrange<BitIter>{bitIter(0), bitIter(derived().size())};
    }


    // ** default implementations for rank and select **
    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rank(Index pos) const noexcept {
        if constexpr (IsOne) {
            return derived().rankOne(pos);
        } else {
            return derived().rankZero(pos);
        }
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankZero(Index pos) const { return pos - derived().rankOne(pos); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankZeroUnchecked(Index pos) const noexcept {
        return pos - derived().rankOneUnchecked(pos);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOne(Index pos) const {
        if (pos >= derived().numBits() || pos < 0) [[unlikely]] {
            throw std::invalid_argument("invalid position for rank query");
        }
        return derived().rankOneUnchecked(pos);
    }

    // This implementation is very slow and should be overwritten in a derived class
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOneUnchecked(Index pos) const {
        return derived().template rank<true>(pos);
    }

    /// Return the position `i` such that `getBit(i)` returns the `rank`th one, counting from 0.
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectOne(Index rank) const {
        if (rank < 0 || rank >= derived().size()) [[unlikely]] {
            throw std::invalid_argument("invalid rank for select one query: " + std::to_string(rank));
        }
        return derived().selectOneUnchecked(rank);
    }

    /// Return the position `i` such that `getBit(i)` returns the `rank`th zero, counting from 0
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectZero(Index rank) const {
        if (rank < 0 || rank >= derived().size()) [[unlikely]] {
            throw std::invalid_argument("invalid rank for select zero query: " + std::to_string(rank));
        }
        return derived().selectZeroUnchecked(rank);
    }

    /// Return the position `i` such that `getBit(i)` returns the `rank`th one, counting from 0.
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectOneUnchecked(Index rank) const {
        return derived().template select<true>(rank);
    }

    /// Return the position `i` such that `getBit(i)` returns the `rank`th zero, counting from 0
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectZeroUnchecked(Index rank) const {
        return derived().template select<false>(rank);
    }
    // This implementation is very slow and should be overwritten in a derived class.
    template<bool IsOne>
    [[nodiscard]] constexpr Index select(Index bitRank) const {
        if constexpr (IsOne) {
            return derived().selectOne(bitRank);
        } else {
            return derived().selectZero(bitRank);
        }
    }

    [[nodiscard]] constexpr std::pair<Index, Index> selectOneAndPrevOne(Index rankOfSecond) const noexcept {
        ADS_ASSUME(rankOfSecond >= 0);
        ADS_ASSUME(rankOfSecond < derived().numOnes());
        if (rankOfSecond == 0) [[unlikely]] {
            return {0, derived().selectOne(rankOfSecond)};
        }
        // for arbitrary bitvectors, there is no better implementation
        return {derived().selectOne(rankOfSecond - 1), derived().selectOne(rankOfSecond)};
    };

    // ** iterators and views over rank and select. Note that no effort is made to speed up consecutive rank calls
    // etc, this simply call rank/select with the current integer value **

    template<bool IsOne>
    using RankIter = RandAccessIter<Derived, decltype(rankFunc<IsOne>)>;

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR RankIter<IsOne> rankIter(Index i) const {
        return RankIter<IsOne>(derived(), rankFunc<IsOne>, i);
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<RankIter<IsOne>> rankView() const noexcept {
        return Subrange<RankIter<IsOne>>{
                derived().template rankIter<IsOne>(0), derived().template rankIter<IsOne>(derived().size())};
    }

    template<bool IsOne>
    using SelectIter = RandAccessIter<Derived, decltype(selectFunc<IsOne>)>;

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR SelectIter<IsOne> selectIter(Index i) const {
        return SelectIter<IsOne>(derived(), selectFunc<IsOne>, i);
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<SelectIter<IsOne>> selectView() const noexcept {
        return Subrange<SelectIter<IsOne>>{derived().template selectIter<IsOne>(0),
                derived().template selectIter<IsOne>(derived().template numBitsEqualTo<IsOne>())};
    }

    // ** Fallback implementations of rank and select that can be used for testing or bitvectors that don't store
    // one type of metadata. They simply do a binary search over select/rank, which is usually not the best strategy. **
    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index inefficientRank(Index i) const {
        auto selects = derived().template selectView<IsOne>();
        return std::lower_bound(selects.begin(), selects.end(), i).index();
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index inefficientSelect(Index bitRank) const {
        auto ranks = derived().template rankView<IsOne>();
        return std::lower_bound(ranks.begin(), ranks.end(), bitRank + 1).index() - 1;
    }

private:
    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankFallback([[maybe_unused]] Index bitRank) const {
#ifdef ADS_NO_FALLBACKS // ensure at compile time that slow operations aren't called
        // static_assert condition must depend on IsOne, else it will fail even if the template isn't being instantiated
        static_assert(IsOne != IsOne, "rank fallback isn't enabled");
        return -1;
#else
        return this->template inefficientRank<IsOne>(bitRank);
#endif
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectFallback([[maybe_unused]] Index bitRank) const {
#ifdef ADS_NO_FALLBACKS // ensure at compile time that slow operations aren't called
        static_assert(IsOne != IsOne, "select fallback isn't enabled");
        return -1;
#else
        return this->template inefficientSelect<IsOne>(bitRank);
#endif
    }

public:
    // ** Printing the bitvector. Again, not the most efficient implementation possible **

    friend std::ostream& operator<<(std::ostream& os, const BitvecBase& bv) {
        std::copy(bv.bitView().begin(), bv.bitView().end(), std::ostream_iterator<bool>(os, ""));
        return os;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR std::string toString() const noexcept {
        std::string res;
        for (bool b : derived().bitView()) {
            if (b) {
                res += '1';
            } else {
                res += '0';
            }
        }
        return res;
    }
};

} // namespace ads


#endif // ADS_BITVEC_HPP
