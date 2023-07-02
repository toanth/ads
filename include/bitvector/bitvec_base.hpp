#ifndef ADS_BITVEC_BASE_HPP
#define ADS_BITVEC_BASE_HPP

#include "../bit.hpp"
#include "../bit_access.hpp"
#include "../common.hpp"
#include "../rand_access_iter.hpp"

namespace ads {

#ifdef ADS_HAS_CPP20
/// Although all the following methods must be present, inheriting from BitvecBase provides default implementations for most.
template<typename T>
concept IsBitvec = requires(T& t, const T& ct) {
    T();
    T(Index(), Limb());
    T(Index(), Limb(), (Limb*)nullptr);
    T("");
    T("", Index());
    T("", Index(), (Limb*)nullptr);
    typename T::BlockRank;
    typename T::SuperblockRank;
    { T::superblockSize() } -> std::convertible_to<Index>;
    { T::blockSize() } -> std::convertible_to<Index>;
    { T::numLimbsInSuperblock() } -> std::convertible_to<Index>;
    { T::numLimbsInBlock() } -> std::convertible_to<Index>;
    { T::bytesPerSuperblockRank() } -> std::convertible_to<Index>;
    { T::bytesPerBlockRank() } -> std::convertible_to<Index>;
    { T::numBlocksInSuperblock() } -> std::convertible_to<Index>;
    { T::numBlocksForBits(Index()) } -> std::convertible_to<Index>;
    { T::numSuperblocksForBits(Index()) } -> std::convertible_to<Index>;
    { T::allocatedSizeInLimbsForLimbs(Index()) } -> std::convertible_to<Index>;
    { T::allocatedSizeInLimbsForBits(Index()) } -> std::convertible_to<Index>;
    { ct.allocatedSizeInLimbs() } -> std::convertible_to<Index>;
    { ct.sizeInLimbs() } -> std::convertible_to<Index>;
    { ct.sizeInBits() } -> std::convertible_to<Index>;
    { ct.size() } -> std::convertible_to<Index>;
    { ct.getLimb(Index()) } -> std::convertible_to<const U64&>;
    { ct.getSuperblockRank(Index()) } -> std::convertible_to<const U64&>;
    t.setSuperblockRank(Index(), Index());
    { ct.getBlockRank(Index()) } -> std::convertible_to<Index>;
    t.setBlockRank(Index(), Index());
    t.setBlockRank_(Index(), Index(), Index());
    { ct.numLimbs() } -> std::convertible_to<Index>;
    { ct.numBlocks() } -> std::convertible_to<Index>;
    { ct.numSuperblocks() } -> std::convertible_to<Index>;
};
#define ADS_BITVEC_CONCEPT IsBitvec
#else
#define ADS_BITVEC_CONCEPT class
#endif



/// \brief CRTP base class of all bitvectors.
/// This works like "normal" inheritance with virtual methods, except that "virtual dispatch" 1. is opt-in by calling
/// derived().func(), 2. works for static methods and templates, 3. is resolved at compile time, and 4. allows
/// signatures to differ. On the downside, there is no single baseclass, so there is no built-in type erasure over
/// different bitvector implementations. The main benefit over simply using a concept is that this allows default
/// implementations for all methods. The `derived()` calls are admittedly ugly but would be unnecessary with C++23's
/// deduced this.
/// \tparam Derived The (most) derived Bitvector class, which needs to inherit from this template like
/// this: `class ActualBitvector : BitvecBase<ActualBitvector>`. Inheritance can be private; the base should be declared
/// as a friend of the derived class.
template<typename Derived, typename BlockRankImpl, typename SuperblockRankImpl = U64> // Don't require Derived to model IsBitvec to avoid recursively depending on this implementation
class [[nodiscard]] BitvecBase {
    friend Derived;

    // ** These aliases can be used unqualified (ie without Derived::) because they should refer to the actual types **
    using BlockRank = BlockRankImpl;
    using SuperblockRank = SuperblockRankImpl;

    /// This owns the uninitialized allocated memory. Being in the base class ensures the memory is (de)allocated at the right time.
    /// In many cases, such as for the Elias-Fano class, this will hold nullptr because the actual owner is another object.
    Allocation<> allocation = Allocation();

    // ** The derived() methods allow simulation of virtual dispatch at compile time **
    // Inlining these functions even in debug builds improves speed and debugging experience
    ADS_FORCE_INLINE([[nodiscard]] ADS_CPP20_CONSTEXPR Derived& derived() noexcept) {
        return static_cast<Derived&>(*this);
    }
    ADS_FORCE_INLINE([[nodiscard]] ADS_CPP20_CONSTEXPR const Derived& derived() const noexcept) {
        return static_cast<const Derived&>(*this);
    }

    /// \brief Allocates memory for the Bitvector or uses the provided ptr, assuming it points to enough memory.
    /// \param numBits The number of bits stored in this bitvector, usually less than the total allocated space.
    /// \param ptr If nullptr, allocate memory. Else, use this pointer; it must point to enough memory.
    ADS_CPP20_CONSTEXPR BitvecBase(Index numBits, Limb* ptr = nullptr) noexcept
        : allocation(Derived::allocatedSizeInLimbsForBits(numBits), ptr) {}

public:
    constexpr BitvecBase() noexcept = default;

    static ADS_CPP20_CONSTEXPR Derived uninitializedForSize(Index numBits, Limb* mem = nullptr) noexcept {
        return Derived(UninitializedTag{}, numBits, mem);
    }

    ADS_CPP20_CONSTEXPR void fill(Limb value) noexcept {
        for (Index i = 0; i + 1 < derived().sizeInLimbs(); ++i) {
            derived().setLimb(i, value);
        }
        Index shift = 64 * derived().sizeInLimbs() - derived().sizeInBits();
        ADS_ASSUME(shift >= 0);
        ADS_ASSUME(shift < 64);
        value = (value << shift) >> shift;
        derived().setLimb(derived().sizeInLimbs() - 1, value);
        derived().buildMetadata();
    }

    ADS_CPP20_CONSTEXPR void copyFrom(Span<const Limb> values) noexcept {
        ADS_ASSUME(values.size() == sizeInLimbs());
        for (Index i = 0; i + 1 < derived().sizeInLimbs(); ++i) {
            derived().setLimb(i, values[i]);
        }
        Index shift = 64 * derived().sizeInLimbs() - derived().sizeInBits();
        ADS_ASSUME(shift >= 0);
        ADS_ASSUME(shift < 64);
        Limb value = (values[values.size() - 1] << shift) >> shift;
        derived().setLimb(derived().sizeInLimbs() - 1, value);
        derived().buildMetadata();
    }

private:
    // *** Default implementations for private Bitvector operations. ***

    // ** Operations useful for implementing getters and setters, not called directly in this class.
    // Note that a Derived class may not need all of these and therefore may not implement them, in which case
    // it should overwrite the getters and setters instead. **

    [[nodiscard]] static ADS_CPP20_CONSTEXPR Index limbIdxInArray(Index i) noexcept { return i; }

    [[nodiscard]] static ADS_CPP20_CONSTEXPR Index blockRankIdxInArray(Index i) noexcept { return i; }

    [[nodiscard]] static ADS_CPP20_CONSTEXPR Index superblockRankIdxInArray(Index i) noexcept { return i; }

    // ** Not implementing either these operations or the getters/setters will cause a stack overflow due to infinite recursion **
    [[nodiscard]] ADS_CPP20_CONSTEXPR const Limb* getLimbArray() const noexcept { return derived().getLimbArray(); }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb* getLimbArray() noexcept { return derived().getLimbArray(); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const View<BlockRank>& getBlockRankArray() const noexcept {
        return derived().getBlockRankArray();
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR View<BlockRank>& getBlockRankArray() noexcept {
        return derived().getBlockRankArray();
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const View<SuperblockRank>& getSuperblockRankArray() const noexcept {
        return derived().getSuperblockRankArray();
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR View<SuperblockRank>& getSuperblockRankArray() noexcept {
        return derived().getSuperblockRankArray();
    }

    // ** An implementation of this interface could ignore the previous methods and only override the following **

    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb& getLimbRef(Index i) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().sizeInLimbs());
        Index idx = derived().limbIdxInArray(i);
        ADS_ASSUME(idx >= 0);
        ADS_ASSUME(idx < allocation.size());
        return derived().getLimbArray()[idx];
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const Limb& getLimbRef(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().sizeInLimbs());
        Index idx = derived().limbIdxInArray(i);
        ADS_ASSUME(idx >= 0);
        ADS_ASSUME(idx < allocation.size());
        return derived().getLimbArray()[idx];
    }

public:
    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb getLimb(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().sizeInLimbs());
        return derived().getLimbRef(i);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR BlockRank getBlockRank(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numBlocks());
        return derived().getBlockRankArray()[derived().blockRankIdxInArray(i)];
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR SuperblockRank getSuperblockRank(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        // the number of superblock ranks is one larger than the number of superblocks
        ADS_ASSUME(i <= derived().numSuperblocks());
        return derived().getSuperblockRankArray()[derived().superblockRankIdxInArray(i)];
    }


    ADS_CPP20_CONSTEXPR void setLimb(Index i, Limb newVal) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().sizeInLimbs());
        derived().getLimbRef(i) = newVal;
    }

    ADS_CPP20_CONSTEXPR void setBlockRank(Index i, BlockRank newVal) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numBlocks());
        derived().getBlockRankArray().setBits(derived().blockRankIdxInArray(i), newVal);
    }

    ADS_CPP20_CONSTEXPR void setBlockRank_(Index superblockIdx, Index blockInSuperblock, BlockRank newVal) noexcept {
        ADS_ASSUME(superblockIdx >= 0);
        ADS_ASSUME(superblockIdx < derived().numSuperblocks());
        ADS_ASSUME(blockInSuperblock >= 0);
        ADS_ASSUME(blockInSuperblock < derived().numBlocksInSuperblock());
        derived().setBlockRank(superblockIdx * derived().numBlocksInSuperblock() + blockInSuperblock, newVal);
    }
    ADS_CPP20_CONSTEXPR void setSuperblockRank(Index i, SuperblockRank newVal) noexcept {
        ADS_ASSUME(i >= 0);
        // the number of superblock ranks is one larger than the number of superblocks
        ADS_ASSUME(i <= derived().numSuperblocks());
        derived().getSuperblockRankArray().setBits(derived().superblockRankIdxInArray(i), newVal);
    }

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

    constexpr static auto getBitFunc = [](const auto& bv, Index i) -> bool { return bv.getBit(i); };
    constexpr static auto getLimbFunc = [](const auto& bv, Index i) -> Limb { return bv.getLimb(i); };
    template<bool IsOne>
    constexpr static auto rankFunc = [](const auto& bv, Index i) -> bool { return bv.template rank<IsOne>(i); };
    template<bool IsOne>
    constexpr static auto selectFunc = [](const auto& bv, Index i) -> bool { return bv.template select<IsOne>(i); };

    ADS_CPP20_CONSTEXPR void initFromStr(std::string_view str, Index base = 2) {
        assert(derived().sizeInBits() == str.size() * intLog2(base)); // should only be called after allocation
        if (base != 2 && base != 4 && base != 16) [[unlikely]] {
            throw std::invalid_argument("base must be one of 2, 4, or 16");
        }
        const Index log2ofBase = intLog2(base);
        const std::size_t charsPerLimb = 64 / log2ofBase;
        Index i = 0;
        for (Index superblock = 0; superblock < derived().numSuperblocks(); ++superblock) {
            for (Index j = 0; j < derived().numLimbsInSuperblock(); ++j) {
                if (str.empty()) break;
                std::size_t numToParse = std::min(charsPerLimb, str.size());
                std::string_view toParse = str.substr(0, numToParse);
                std::uint64_t res;
                auto err = fromChars(toParse.data(), toParse.data() + toParse.size(), res, int(base));
                if (err.ec != std::errc()) [[unlikely]] {
                    throw std::invalid_argument(std::make_error_code(err.ec).message());
                } else if (err.ptr != toParse.data() + toParse.size()) [[unlikely]] {
                    throw std::invalid_argument("invalid character found");
                }
                if (numToParse < charsPerLimb) {
                    res <<= (charsPerLimb - numToParse) * log2ofBase;
                }
                str.remove_prefix(numToParse);
                res = reverseBits(res);
                derived().setLimb(i, res);
                ++i;
            }
            derived().buildRankMetadata(superblock);
        }
        derived().finalizeRankMetadata();
        derived().buildSelectMetadata();
        assert(str.empty());
    }

public:
    // *** Default implementations for (static) getters ***

    // ** Getting the size in bits or limbs. Usually, it's enough to implement size() in the derived class. **

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return derived().sizeInBits(); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index sizeInBits() const noexcept { return derived().size(); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index sizeInLimbs() const noexcept { return derived().numLimbs(); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numLimbs() const noexcept { return roundUpDiv(derived().size(), 64); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBlocks() const noexcept {
        return Derived::numBlocksForBits(derived().size());
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numSuperblocks() const noexcept {
        return Derived::numSuperblocksForBits(derived().size());
    }

    [[nodiscard]] static constexpr Index superblockSize() noexcept { return Derived::numLimbsInSuperblock() * 64; }

    [[nodiscard]] static constexpr Index numLimbsInSuperblock() noexcept {
        return roundUpDiv(Derived::superblockSize(), 64);
    }

    [[nodiscard]] static constexpr Index numLimbsInBlock() noexcept { return roundUpDiv(Derived::blockSize(), 64); }

    [[nodiscard]] static constexpr Index blockSize() noexcept { return Derived::numLimbsInBlock() * 64; }

    [[nodiscard]] static constexpr Index numBlocksInSuperblock() noexcept {
        return Derived::numLimbsInSuperblock() / Derived::numLimbsInBlock();
    }

    [[nodiscard]] static constexpr Index bytesPerSuperblockRank() noexcept {
        return roundUpLog2(Derived::superblockSize());
    }

    [[nodiscard]] static constexpr Index bytesPerBlockRank() noexcept { return roundUpLog2(Derived::blockSize()); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index allocatedSizeInLimbs() const noexcept {
        return Derived::allocatedSizeInLimbsForLimbs(derived().sizeInLimbs());
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index allocatedSizeInBits() const noexcept {
        return derived().allocatedSizeInLimbs() * 64; // only complete limbs are allocated
    }

    [[nodiscard]] static constexpr Index allocatedSizeInLimbsForBits(Index numBits) noexcept {
        return Derived::allocatedSizeInLimbsForLimbs(roundUpDiv(numBits, 64));
    }

    [[nodiscard]] static constexpr Index allocatedSizeInLimbsForLimbs(Index numLimbs) noexcept {
        return Derived::allocatedSizeInLimbsForBits(numLimbs * 64);
    }

    [[nodiscard]] static constexpr Index numBlocksForBits(Index numBits) noexcept {
        return roundUpDiv(numBits, Derived::blockSize());
    }

    [[nodiscard]] static constexpr Index numSuperblocksForBits(Index numBits) noexcept {
        return roundUpDiv(numBits, Derived::superblockSize());
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numOnes() const noexcept {
        return derived().getSuperblockRank(derived().numSuperblocks());
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numZeros() const noexcept {
        return derived().sizeInBits() - derived().numOnes();
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numBitsEqualTo() const noexcept {
        if constexpr (IsOne) {
            return derived().numOnes();
        } else {
            return derived().numZeros();
        }
    }


    // ** The following functions provide a view over bits or limbs. The bitView isn't the most efficient **
    using BitIter = RandAccessIter<Derived, decltype(getBitFunc)>;

    [[nodiscard]] ADS_CPP20_CONSTEXPR BitIter bitIter(Index i) const { return BitIter(derived(), getBitFunc, i); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<BitIter> bitView() const noexcept {
        return Subrange<BitIter>{bitIter(0), bitIter(derived().size())};
    }

    using LimbIter = RandAccessIter<Derived, decltype(getLimbFunc)>;

    [[nodiscard]] ADS_CPP20_CONSTEXPR LimbIter limbIter(Index i) const { return LimbIter(derived(), getLimbFunc, i); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<LimbIter> limbView() const noexcept {
        return Subrange<LimbIter>{limbIter(0), limbIter(derived().sizeInLimbs())};
    }

    // ** A derived class may change the return type of this function, eg to Subrange<Limb>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<Limb> superblockLimbs(Index superblockIdx) noexcept {
        Index size = superblockIdx == derived().numSuperblocks() - 1 ?
                             (derived().numLimbs() - 1) % derived().numLimbsInSuperblock() + 1 :
                             derived().numLimbsInSuperblock();
        return Span<Limb>(&derived().getLimbRef(superblockIdx * derived().numLimbsInSuperblock()), size);
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<const Limb> superblockLimbs(Index superblockIdx) const noexcept {
        Index size = superblockIdx == derived().numSuperblocks() - 1 ?
                             (derived().numLimbs() - 1) % derived().numLimbsInSuperblock() + 1 :
                             derived().numLimbsInSuperblock();
        return Span<const Limb>(&derived().getLimbRef(superblockIdx * derived().numLimbsInSuperblock()), size);
    }


    // *** Default implementations for Bitvector operations that are more than just (static) getters/setters ***


    // ** low-level operations that are nevertheless public because sometimes it's convenient or
    // faster to let the user call them **

    ADS_CPP20_CONSTEXPR void buildRankMetadata(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < derived().numSuperblocks());
        if (superblockIdx == 0) {
            derived().setSuperblockRank(0, 0);
        }
        auto s = derived().superblockLimbs(superblockIdx);
        Index inSuperblockSoFar = 0;
        for (Index i = 0; i < s.size(); ++i) {
            if (i % derived().numLimbsInBlock() == 0) {
                //    std::cout << "setting block count for block " << i << " to " << inSuperblockSoFar << ", s[i] is "
                //              << std::hex << s[i] << std::dec << ", popcount " << popcount(s[i]) << std::endl;
                derived().setBlockRank_(superblockIdx, i / derived().numLimbsInBlock(), inSuperblockSoFar);
            }
            inSuperblockSoFar += popcount(s[i]);
        }
        assert(inSuperblockSoFar <= derived().superblockSize());
        if (superblockIdx + 1 < derived().numSuperblocks()) {
            assert(s.size() % derived().numLimbsInBlock() == 0);
        }
        // there are numSuperblocks() + 1 super block counts
        derived().setSuperblockRank(superblockIdx + 1, derived().getSuperblockRank(superblockIdx) + inSuperblockSoFar);
    }

    ADS_CPP20_CONSTEXPR void finalizeRankMetadata() noexcept {
        // do nothing
    }

    ADS_CPP20_CONSTEXPR void buildRankMetadata() noexcept {
        for (Index i = 0; i < derived().numSuperblocks(); ++i) {
            derived().buildRankMetadata(i);
        }
        derived().finalizeRankMetadata();
    }

    ADS_CPP20_CONSTEXPR void buildSelectMetadata() noexcept {
        // nothing to do in this default implementation
    }

    ADS_CPP20_CONSTEXPR void buildMetadata() noexcept {
        derived().buildRankMetadata();
        derived().buildSelectMetadata();
    }

    // ** default implementations for rank and select **

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankZero(Index pos) const { return pos - derived().rankOne(pos); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankZeroUnchecked(Index pos) const noexcept {
        return pos - derived().rankOneUnchecked(pos);
    }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOne(Index pos) const {
        if (pos >= derived().sizeInBits() || pos < 0) [[unlikely]] {
            throw std::invalid_argument("invalid position for rank query");
        }
        return derived().rankOneUnchecked(pos);
    }


    /// Unlike rankOne(), this dDoesn't check that `pos` is valid, although that gives close to no measurable performance benefits.
    /// However, the combined ASSUME macros do improve performance by quite a bit (if the compiler couldn't assume that pos >= 0,
    //// performance would actually be significantly lower than with the throwing checks)
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOneUnchecked(Index pos) const noexcept {
        ADS_ASSUME(0 <= pos);
        ADS_ASSUME(pos < derived().sizeInBits());
        if (pos <= 0) [[unlikely]] {
            return 0;
        }
        --pos;
        ADS_ASSUME(pos >= 0);
        Index elemIdx = pos / 64;
        Limb mask = Limb(-1) >> (63 - pos % 64);
        Index superblockIdx = elemIdx / derived().numLimbsInSuperblock();
        ADS_ASSUME(superblockIdx >= 0);
        Index blockIdx = elemIdx / derived().numLimbsInBlock();
        ADS_ASSUME(blockIdx >= 0);
        Index res = derived().getSuperblockRank(superblockIdx) + derived().getBlockRank(blockIdx);
        ADS_ASSUME(res >= 0);
        for (Index i = blockIdx * derived().numLimbsInBlock(); i < elemIdx; ++i) {
            ADS_ASSUME(elemIdx - i < derived().numLimbsInBlock());
            res += popcount(derived().getLimb(i));
        }
        return res + popcount(derived().getLimb(elemIdx) & mask);
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rank(Index pos) const noexcept {
        if constexpr (IsOne) {
            return derived().rankOne(pos);
        } else {
            return derived().rankZero(pos);
        }
    }


    /// Return the position `i` such that `getBit(i)` returns the `rank`th one, counting from 0.
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectOne(Index rank) const {
        return derived().template select<true>(rank);
    }

    /// Return the position `i` such that `getBit(i)` returns the `rank`th zero, counting from 0
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectZero(Index rank) const {
        return derived().template select<false>(rank);
    }

private:
    // ** internal implementations for the default select. They simply do binary searches, but the
    // constant factor is smaller than for inefficientSelect(). Still, some Bitvector implementations
    // provide faster select operations ** // TODO: Actuallu implement
    template<bool IsOne>
    [[nodiscard]] constexpr Index selectSuperBlockIdx(Index& bitRank) const {
        constexpr Index linearFallbackSize = 8;
        auto rankFunc = [this](Index i) noexcept {
            ADS_ASSUME(i >= 0);
            ADS_ASSUME(i <= derived().numSuperblocks());
            Index rankOne = derived().getSuperblockRank(i);
            ADS_ASSUME(rankOne >= 0);
            if constexpr (IsOne) {
                return rankOne;
            } else {
                // if i is derived().numSuperblocks(), this can be greater than the number of zeros in the bitvector,
                // but that's not a problem
                Index numBitsBefore = i * derived().superblockSize();
                ADS_ASSUME(rankOne <= numBitsBefore);
                return numBitsBefore - rankOne;
            }
        };
        ADS_ASSUME(bitRank >= 0);
        // + 1 because we're searching for the first superblock where the rank is greater than bitRank
        Index l = bitRank / derived().superblockSize() + 1;
        Index u = derived().numSuperblocks();
        ADS_ASSUME(l > 0);
        ADS_ASSUME(l <= u);
        if (u - l > linearFallbackSize) {
            // set u close to the expected location for iid ones with 50% probability, then increase exponentially
            // until it is an upper bound. Unlike binary search, this starts with a less pessimistic search window and
            // should hopefully be easier on the branch predictor. This improves performance for random values but hurts for especially hard cases.
            u = l;
            do {
                l = u;
                u *= 2;
                if (u >= derived().numSuperblocks()) {
                    u = derived().numSuperblocks();
                    break;
                }
            } while (rankFunc(u) <= bitRank);
            while (u - l > linearFallbackSize) {
                ADS_ASSUME(0 < l);
                ADS_ASSUME(l <= u);
                Index mid = (l + u) / 2;
                Index midRank = rankFunc(mid);
                ADS_ASSUME(midRank >= rankFunc(l));
                if (midRank <= bitRank) {
                    l = mid;
                } else {
                    u = mid;
                }
            }
        }
        ADS_ASSUME(u - l <= linearFallbackSize);
        ADS_ASSUME(l <= u);
        for (Index i = l; i < u; ++i) {
            if (rankFunc(i) > bitRank) {
                bitRank -= rankFunc(i - 1);
                return i - 1;
            }
        }
        ADS_ASSUME(rankFunc(u) >= bitRank);
        bitRank -= rankFunc(u - 1);
        return u - 1;
    }

    // TODO: For most bitvectors, superblock indices can be represented with 32 bit values.
    // TODO: Also, it may be worth investigating if binary search works better if the range isn't split in two equal
    // sized halves but instead the requested bitRank as a fraction of the total rank in the (bv|superblock|block) is
    // used to generate the pivot element. The problem of this approach is that it can be very inefficient for
    // adversarial patters such as a superblock with 100 ones, all of which are in the last 100 bits and a selectOne(0)
    // query. Maybe some combination with splitting in equal-sized halves can give reasonable worst-case guarantees?
    template<bool IsOne>
    [[nodiscard]] constexpr Index selectBlockIdx(Index& bitRank, Index superBlockIdx) const noexcept {
        auto rankFunc = [this](Index i) noexcept {
            Index rankOne = derived().getBlockRank(i);
            if constexpr (IsOne) {
                return rankOne;
            } else {
                Index numBitsBefore = (i % derived().numBlocksInSuperblock()) * derived().blockSize();
                ADS_ASSUME(rankOne <= numBitsBefore);
                return numBitsBefore - rankOne;
            }
        };
        constexpr Index linearFallbackSize = 2 * sizeof(Limb) / sizeof(BlockRank);
        ADS_ASSUME(superBlockIdx >= 0);
        ADS_ASSUME(superBlockIdx < derived().numSuperblocks());
        ADS_ASSUME(bitRank >= 0);
        // we're searching for the first block with count strictly greater than bitRank
        Index l = superBlockIdx * derived().numBlocksInSuperblock() + 1;
        Index u = std::min(l + derived().numBlocksInSuperblock() - 1, derived().numBlocks());
        ADS_ASSUME(u >= l);
        ADS_ASSUME(u - l < derived().numBlocksInSuperblock());
        while (u - l > linearFallbackSize) {
            Index mid = (l + u) / 2;
            Index midRank = rankFunc(mid);
            if (midRank > bitRank) {
                u = mid;
            } else {
                l = mid;
            }
        }
        ADS_ASSUME(u >= l);
        ADS_ASSUME(u - l <= linearFallbackSize);
        ADS_ASSUME(l > superBlockIdx * derived().numBlocksInSuperblock());
        for (Index i = l; i < u; ++i) {
            ADS_ASSUME(i == 0 || derived().getBlockRank(i) >= derived().getBlockRank(i - 1));
            if (rankFunc(i) > bitRank) {
                bitRank -= rankFunc(i - 1);
                return i - 1;
            }
        }
        bitRank -= rankFunc(u - 1);
        return u - 1;
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectLimbIdx(Index& bitRank, Index blockIdx) const noexcept {
        Index first = blockIdx * derived().numLimbsInBlock();
        auto rankFunc = [this](Index i) noexcept {
            Limb l = derived().getLimb(i);
            if constexpr (IsOne) {
                return popcount(l);
            } else {
                return 64 - popcount(l);
            }
        };
        for (Index i = first; i < numLimbs(); ++i) {
            Index rank = rankFunc(i);
            ADS_ASSUME(i - first < derived().numLimbsInBlock());
            ADS_ASSUME(rank >= 0);
            if (rank > bitRank) {
                return i;
            }
            bitRank -= rank;
            ADS_ASSUME(bitRank >= 0);
        }
        return numLimbs() - 1;
    }

    template<bool IsOne>
    [[nodiscard]] constexpr Index selectBitIdx(Limb limb, Index bitIndex) const noexcept {
        if constexpr (IsOne) {
            return u64Select(limb, bitIndex);
        } else {
            return u64Select(~limb, bitIndex);
        }
    }

public:
    template<bool IsOne>
    [[nodiscard]] constexpr Index select(Index bitRank) const {
        if (bitRank < 0 || bitRank >= derived().size()) [[unlikely]] {
            throw std::invalid_argument("invalid rank for select query: " + std::to_string(bitRank));
        }
        Index superBlockIdx = derived().template selectSuperBlockIdx<IsOne>(bitRank);
        Index blockIdx = derived().template selectBlockIdx<IsOne>(bitRank, superBlockIdx);
        Index limbIdx = derived().template selectLimbIdx<IsOne>(bitRank, blockIdx);
        Index bitIdx = derived().template selectBitIdx<IsOne>(derived().getLimb(limbIdx), bitRank);
        return limbIdx * 64 + bitIdx;
    }


    // ** iterators and views over rank and select. Note that no effort is made to speed up consecutive rank calls etc,
    // this simply call rank/select with the current integer value **

    template<bool IsOne>
    using RankIter = RandAccessIter<Derived, decltype(rankFunc<IsOne>)>;

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR LimbIter rankIter(Index i) const {
        return RankIter<IsOne>(derived(), rankFunc<IsOne>, i);
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<RankIter<IsOne>> rankView() const noexcept {
        return Subrange<RankIter<IsOne>>{rankIter<IsOne>(0), rankIter<IsOne>(derived().size())};
    }

    template<bool IsOne>
    using SelectIter = RandAccessIter<Derived, decltype(selectFunc<IsOne>)>;

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR LimbIter selectIter(Index i) const {
        return SelectIter<IsOne>(derived(), selectFunc<IsOne>, i);
    }

    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Subrange<SelectIter<IsOne>> selectView() const noexcept {
        return Subrange<SelectIter<IsOne>>{selectIter<IsOne>(0), selectIter<IsOne>(derived().template numBitsEqualTo<IsOne>())};
    }

    // ** printing the bitvector **

    friend std::ostream& operator<<(std::ostream& os, const Derived& bv) {
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

    // ** Fallback implementations of rank and select that can be used for testing or bitvectors that don't store
    // metadata. They simply do a binary search over select/rank, which is usually not the best strategy. **
    template<bool IsOne>
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index inefficientRank(Index i) const {
        if (i < 0 || i >= derived().sizeInBits()) [[unlikely]] {
            throw std::invalid_argument("invalid position for rank query: " + std::to_string(i));
        }
        auto selects = derived().template selectView<IsOne>();
        return std::lower_bound(selects.begin(), selects.end(), i).index();
    }

    template<bool IsOne> // TODO: Test
    [[nodiscard]] ADS_CPP20_CONSTEXPR Index inefficientSelect(Index bitRank) const {
        if (bitRank < 0 || bitRank >= derived().template numBitsEqualTo<IsOne>()) [[unlikely]] {
            throw std::invalid_argument("invalid rank for select query: " + std::to_string(bitRank));
        }
        auto ranks = derived().template rankView<IsOne>();
        return std::lower_bound(ranks.begin(), ranks.end(), bitRank).index();
    }
};


template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] ADS_CPP20_CONSTEXPR bool operator==(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return lhs.sizeInBits() == rhs.sizeInBits()
           && std::equal(lhs.limbView().begin(), lhs.limbView().end(), rhs.limbView().begin());
}
template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] ADS_CPP20_CONSTEXPR bool operator!=(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return !(rhs == lhs);
}

#ifdef ADS_HAS_CPP20

template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] constexpr std::strong_ordering operator<=>(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return lexicographicalCompareThreeWay(
            lhs.limbView().begin(), lhs.limbView().end(), rhs.limbView().begin(), rhs.limbView().end());
}

#else

template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] bool operator<(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return std::lexicographical_compare(
            lhs.limbView().begin(), lhs.limbView().end(), rhs.limbView().begin(), rhs.limbView().end());
}
template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] bool operator>(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return rhs < lhs;
}
template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] bool operator<=(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return !(rhs < lhs);
}
template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] bool operator>=(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return !(lhs < rhs);
}

#endif



} // namespace ads

#endif // ADS_BITVEC_BASE_HPP
