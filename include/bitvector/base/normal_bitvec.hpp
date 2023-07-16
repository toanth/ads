#ifndef ADS_NORMAL_BITVEC_HPP
#define ADS_NORMAL_BITVEC_HPP

#include "bitvec.hpp"

namespace ads {


/// \brief CRTP base class of all normal bitvectors, which are bitvectors that store their bit sequence as an actual
/// sequence of bits, organized in limbs. A fixed number of limbs form a block with implementation-defined meaning. See
/// BitvecBase for a CRTP discussion. Derived classes should declare both this class and its base class, available under
/// the Base alias, as friends. \tparam Derived The actual bitvector, which inherits from NormalBitvecBase<Derived>.
// Don't require Derived to model IsBitvec to avoid recursively depending on this implementation
template<typename Derived, Index NumLimbsInBlock, Index NumLimbsInCacheLine = U64_PER_CACHELINE>
class [[nodiscard]] NormalBitvecBase : public BitvecBase<Derived> {
protected:
    friend Derived;

    static_assert(NumLimbsInCacheLine % NumLimbsInBlock == 0);

    using Base = BitvecBase<Derived>;
    using Base::allocation;
    using Base::Base;
    using Base::derived;

public:
    ADS_CPP20_CONSTEXPR void fill(Limb value) noexcept { copyFrom(repeatView(value, derived().numLimbs())); }

    template<typename LimbRange>
    ADS_CPP20_CONSTEXPR void copyFrom(const LimbRange& values) noexcept {
        ADS_ASSUME(values.size() == derived().numLimbs());
        for (Index i = 0; i + 1 < derived().numLimbs(); ++i) {
            derived().setLimb(i, values.begin()[i]);
        }
        if (derived().numLimbs() > 0) {
            Index shift = 64 * derived().numLimbs() - derived().numBits();
            ADS_ASSUME(shift >= 0);
            ADS_ASSUME(shift < 64);
            Limb value = (values[values.size() - 1] << shift) >> shift;
            derived().setLimb(derived().numLimbs() - 1, value);
        }
        derived().buildMetadata();
    }

protected:
    ADS_CPP20_CONSTEXPR void completeWithZeros() noexcept {
        if (derived().size() % 64 != 0) {
            Limb last = derived().getLimb(derived().numLimbs() - 1);
            Index numBitsToRemove = 64 - (derived().size() % 64);
            last = (last << numBitsToRemove) >> numBitsToRemove;
            derived().setLimb(derived().numLimbs() - 1, last);
        }
        for (Index limb = derived().numLimbs(); limb < derived().numAccessibleLimbs(); ++limb) {
            derived().setLimb(limb, 0);
        }
    }

    ADS_CPP20_CONSTEXPR void initFromStr(std::string_view str, Index base = 2) {
        assert(derived().numBits() == Index(str.size()) * intLog2(base)); // should only be called after allocation
        auto limbValues = this->limbViewFromStringView(str, base);
        assert((limbValues.size() - 1) * 64 < derived().size());
        assert(limbValues.size() * 64 >= derived().size());
        assert(limbValues.size() == derived().numLimbs());
        for (Index i = 0; i < limbValues.size(); ++i) {
            derived().setLimb(i, limbValues[i]);
        }
        derived().buildMetadata();
    }

    // *** Default implementations for private Bitvector operations. ***

    // ** Operations useful for implementing getters and setters, not called directly in this class.
    // Note that a Derived class may not need all of these and therefore may not implement them, in which case
    // it should overwrite the getters and setters instead. **

    [[nodiscard]] static constexpr Index limbIdxInCacheLineArray(Index i) noexcept {
        return i / Derived::numLimbsInCacheLine();
    }
    [[nodiscard]] static constexpr Index limbIdxInCacheLine(Index i) noexcept {
        return i % Derived::numLimbsInCacheLine();
    }

public:
    [[nodiscard]] static constexpr Index numLimbsInCacheLine() noexcept { return NumLimbsInCacheLine; }


    [[nodiscard]] static constexpr Index numBlocksInCacheLine() noexcept {
        Index limbsInCacheLine = Derived::numLimbsInCacheLine();
        Index limbsInBLock = Derived::numLimbsInBlock();
        ADS_ASSUME(limbsInCacheLine % limbsInBLock == 0);
        return limbsInCacheLine / limbsInBLock;
    }

    // ** One of the following two function overload sets must be implemented by the base class **
    //    [[nodiscard]] ADS_CPP20_CONSTEXPR Array<const CacheLine> getCacheLineArray() const noexcept {
    //        return derived().getCacheLineArray();
    //    }
    //    [[nodiscard]] ADS_CPP20_CONSTEXPR Array<CacheLine> getCacheLineArray() noexcept {
    //        return derived().getCacheLineArray();
    //    }

    // or:

    [[nodiscard]] ADS_CPP20_CONSTEXPR const CacheLine& getCompleteCacheLine(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numAccessibleCacheLines());
        return derived().getCacheLineArray()[i];
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR CacheLine& getCompleteCacheLine(Index i) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numAccessibleCacheLines());
        return derived().getCacheLineArray()[i];
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const CacheLine& getCacheLineForLimb(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numAccessibleLimbs());
        i = limbIdxInCacheLineArray(i);
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numAccessibleCacheLines());
        return derived().getCompleteCacheLine(i);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR CacheLine& getCacheLineForLimb(Index i) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numAccessibleLimbs());
        i = limbIdxInCacheLineArray(i);
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numAccessibleCacheLines());
        return derived().getCompleteCacheLine(i);
    }

    // ** An implementation of this interface could ignore the previous methods and only override the following **

    [[nodiscard]] constexpr Span<Limb, NumLimbsInCacheLine> getCacheLine(Index cacheLineIdx) noexcept {
        ADS_ASSUME(cacheLineIdx >= 0);
        ADS_ASSUME(cacheLineIdx < derived().numAccessibleCacheLines());
        CacheLine& cl = derived().getCompleteCacheLine(cacheLineIdx);
        return Span<Limb, NumLimbsInCacheLine>(cl.limbs);
    }

    [[nodiscard]] constexpr Span<const Limb, NumLimbsInCacheLine> getCacheLine(Index cacheLineIdx) const noexcept {
        ADS_ASSUME(cacheLineIdx >= 0);
        ADS_ASSUME(cacheLineIdx < derived().numAccessibleCacheLines());
        const CacheLine& cl = derived().getCompleteCacheLine(cacheLineIdx);
        return Span<const Limb, NumLimbsInCacheLine>(cl.limbs);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb& getLimbRef(Index i) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numAccessibleLimbs());
        Index clIdx = derived().limbIdxInCacheLineArray(i);
        ADS_ASSUME(clIdx >= 0);
        ADS_ASSUME(clIdx < derived().numAccessibleCacheLines());
        Index limbIdx = derived().limbIdxInCacheLine(i);
        ADS_ASSUME(limbIdx >= 0);
        ADS_ASSUME(limbIdx < derived().numLimbsInCacheLine());
        return derived().getCacheLine(clIdx)[limbIdx];
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const Limb& getLimbRef(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numAccessibleLimbs());
        Index clIdx = derived().limbIdxInCacheLineArray(i);
        ADS_ASSUME(clIdx >= 0);
        ADS_ASSUME(clIdx < derived().numAccessibleCacheLines());
        Index limbIdx = derived().limbIdxInCacheLine(i);
        ADS_ASSUME(limbIdx >= 0);
        ADS_ASSUME(limbIdx < derived().numLimbsInCacheLine());
        return derived().getCacheLine(clIdx)[limbIdx];
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb getLimb(Index i) const noexcept { return derived().getLimbRef(i); }

    ADS_CPP20_CONSTEXPR void setLimb(Index i, Limb newVal) noexcept { derived().getLimbRef(i) = newVal; }

    // ** The following two functions rarely need to be overwritten and should be used even less
    // in performance critical code **
    [[nodiscard]] ADS_CPP20_CONSTEXPR bool getBit(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().size());
        return BitwiseAccess<1>::getBits(&derived().getLimbRef(i / 64), i % 64);
    }
    ADS_CPP20_CONSTEXPR void setBit(Index i, bool newVal = true) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().size());
        return BitwiseAccess<1>::setBits(&derived().getLimbRef(i / 64), i % 64, newVal);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Byte getByte(Index idx) const noexcept {
        Limb limb = derived().getLimb(idx / 8);
        return BitwiseAccess<8>::getBits(&limb, 0, idx % 8);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<const Limb, NumLimbsInBlock> getBlock(Index blockIdx) const noexcept {
        Span<const Limb, NumLimbsInCacheLine> cl = getCacheLine(blockIdx / Derived::numBlocksInCacheLine());
        Index offset = (blockIdx % Derived::numBlocksInCacheLine()) * NumLimbsInBlock;
        ADS_ASSUME(offset >= 0);
        ADS_ASSUME(offset < NumLimbsInCacheLine);
        ADS_ASSUME(NumLimbsInCacheLine % NumLimbsInBlock == 0);
        return Span<const Limb, NumLimbsInBlock>(cl.data() + offset);
    }

    template<Group G>
    [[nodiscard]] constexpr auto get(Index idx) const noexcept {
        if constexpr (G == Group::Bit) {
            return derived().getBit(idx);
        } else if constexpr (G == Group::Byte) {
            return derived().getByte(idx);
        } else if constexpr (G == Group::Limb) {
            return derived().getLimb(idx);
        } else if constexpr (G == Group::Block) {
            return derived().getBlock(idx);
        } else if constexpr (G == Group::CacheLine) {
            return derived().getCacheLine(idx);
        } else {
            static_assert(G != G, "Only a SuperblockBitvec uses superblocks");
            return -1;
        }
    }

    template<Group G>
    [[nodiscard]] constexpr Index num() const noexcept {
        return roundUpDiv(derived().size(), Derived::template numBitsIn<G>());
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numCacheLines() const noexcept {
        return derived().template sizeIn<Group::CacheLine>();
    }

    template<Group G>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numAccessible() const noexcept {
        return derived().numCacheLines() * Derived::cacheLineSize() / Derived::template numBitsIn<G>();
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numAccessibleCacheLines() const noexcept {
        return derived().template numAccessible<Group::CacheLine>();
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numAccessibleBlocks() const noexcept {
        return derived().template numAccessible<Group::Block>();
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numLimbs() const noexcept {
        return derived().template sizeIn<Group::Limb>();
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBlocks() const noexcept {
        Index result = derived().template sizeIn<Group::Block>();
        ADS_ASSUME(result == Derived::numBlocksForBits(derived().size()));
        return result;
    }

    /// \brief Unlike numLimbs(), this counts additional limbs that may be appended to make the bitvector size
    /// a multiple of the cache line size or block size. Their values may be indeterminate, depending on the bitvector.
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numAccessibleLimbs() const noexcept { return numAccessible<Group::Limb>(); }

    [[nodiscard]] static constexpr Index numAccessibleBytesForBits(Index numBits) noexcept {
        ADS_ASSUME((Derived::numLimbsInCacheLine() * 64) % Derived::blockSize() == 0);
        return roundUpTo(numBits, CACHELINE_SIZE_BYTES * 8) / 8;
    }

    [[nodiscard]] static constexpr Index numLimbsInBlock() noexcept { return NumLimbsInBlock; }

    [[nodiscard]] static constexpr Index blockSize() noexcept { return Derived::numLimbsInBlock() * 64; }

    [[nodiscard]] static constexpr Index cacheLineSize() noexcept { return Derived::numLimbsInCacheLine() * 64; }

    [[nodiscard]] static constexpr Index bytesPerBlockRank() noexcept { return roundUpLog2(Derived::blockSize()); }


    [[nodiscard]] static constexpr Index numBlocksForBits(Index numBits) noexcept {
        return roundUpDiv(numBits, Derived::blockSize());
    }

    // This function is only necessary when the bitvector is explicitly given memory instead of allocating memory itself
    [[nodiscard]] static constexpr Index requiredAlignment() noexcept { return CACHELINE_SIZE_BYTES; }

    template<Group G>
    using GroupIter = RandAccessIter<Derived, decltype(getGroupFunc<G>)>;

    template<Group G>
    [[nodiscard]] ADS_CPP20_CONSTEXPR GroupIter<G> groupIter(Index i) const {
        return GroupIter<G>(derived(), getGroupFunc<G>, i);
    }

    template<Group G>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<GroupIter<G>> viewOf() const noexcept {
        return Subrange<GroupIter<G>>{groupIter<G>(0), groupIter<G>(derived().template num<G>())};
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<GroupIter<Group::Limb>> limbView() const noexcept {
        return viewOf<Group::Limb>();
    }
};



} // namespace ads


#endif // ADS_NORMAL_BITVEC_HPP
