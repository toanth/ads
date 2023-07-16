#ifndef ADS_TRIVIAL_BITVEC_HPP
#define ADS_TRIVIAL_BITVEC_HPP

#include "../common.hpp"
#include "base/bitvec.hpp"

namespace ads {

enum class Operations { RANK_ONLY, SELECT_ONLY, BOTH };

/// \brief This class precomputes the answers to all possible rank and select queries and stores them in one large array.
/// `getBit(i)` is implemented as `rankOne(i+1) - rankOne(i)` for all i (there are size() + 1 rank entries).
/// For 32 bit rank and select indices, this uses 64 times as much space as the actual bit sequence (which isn't stored),
/// so this class is only useful for small Bitvectors, such as in the second recursion level of the recursive bitvector.
/// For RANK_ONLY or SELECT_ONLY, using 32 bit indices requires 32 times as much space as the actual bit sequence.
template<typename Count = std::uint32_t, Operations Ops = Operations::BOTH>
class TrivialBitvec : public BitvecBase<TrivialBitvec<Count, Ops>> {

    using Base = BitvecBase<TrivialBitvec<Count, Ops>>;
    friend Base;

    Array<Count> ranks = Array<Count>();
    Array<Count> selectAnswers = Array<Count>();

    template<typename LimbRange>
    ADS_CPP20_CONSTEXPR void constructFromLimbRange(const LimbRange& range) noexcept {
        ADS_ASSUME(range.size() * 64 >= size());
        ADS_ASSUME((range.size() - 1) * 64 < size());
        Index previousRank = 0;
        Index numOnesSoFar = 0;
        Index numZerosSoFar = 0;
        Index bitNum = 0;
        for (Limb current : range) {
            for (Index inLimb = 0; inLimb < 64 && bitNum < size(); ++inLimb, ++bitNum) {
                bool bit = BitwiseAccess<1>::getBits(&current, inLimb);
                if constexpr (Ops == Operations::SELECT_ONLY) {
                    ranks.ptr[0] = numOnesSoFar + bit; // ranks has only 1 entry, which counts the total number of ones
                } else {
                    ranks.ptr[bitNum] = previousRank;
                }
                previousRank += bit;
                if constexpr (Ops != Operations::RANK_ONLY) {
                    if (bit) {
                        selectAnswers.ptr[numOnesSoFar++] = bitNum;
                    } else {
                        selectAnswers.ptr[size() - ++numZerosSoFar] = bitNum;
                    }
                }
            }
        }
        if constexpr (Ops != Operations::SELECT_ONLY) {
            ranks.ptr[size()] = previousRank;
        }
        buildMetadata();
    }

    ADS_CPP20_CONSTEXPR static Index numRankEntries(Index numBits) noexcept {
        if constexpr (Ops == Operations::SELECT_ONLY) {
            return 1;
        } else {
            return numBits + 1;
        }
    }

    template<typename Underlying = CacheLine>
    ADS_CPP20_CONSTEXPR TrivialBitvec(UninitializedTag, Index numBits, Underlying* mem) noexcept
        : Base(numBits, mem), ranks(this->allocation.memory(), numRankEntries(numBits)),
          selectAnswers(ranks.ptr + ranks.numT, Ops == Operations::RANK_ONLY ? 0 : numBits) {
        if constexpr (Ops == Operations::SELECT_ONLY) {
            ranks.ptr[0] = 0; // special entry that counts the number of ones
        }
    }

public:
    constexpr TrivialBitvec() noexcept = default;
    ADS_CPP20_CONSTEXPR TrivialBitvec(Index size, Limb fill, CacheLine* mem = nullptr) noexcept
        : TrivialBitvec(UninitializedTag{}, size, mem) {
        constructFromLimbRange(repeatView(fill, roundUpDiv(size, 64)));
    }

    explicit ADS_CPP20_CONSTEXPR TrivialBitvec(Span<const Limb> limbs, CacheLine* mem = nullptr) noexcept
        : TrivialBitvec(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR TrivialBitvec(Span<const Limb> limbs, Index numBits, CacheLine* mem = nullptr) noexcept
        : TrivialBitvec(UninitializedTag{}, numBits, mem) {
        constructFromLimbRange(limbs);
    }


    ADS_CPP20_CONSTEXPR TrivialBitvec(std::string_view str, Index base = 2, CacheLine* mem = nullptr)
        : TrivialBitvec(UninitializedTag{}, str.size() * intLog2(base), mem) {
        constructFromLimbRange(this->limbViewFromStringView(str, base));
    }

    [[nodiscard]] static ADS_CPP20_CONSTEXPR Index allocatedSizeInBytesForBits(Index numBits) noexcept {
        const Index factor = (Ops == Operations::BOTH ? 2 : 1);
        return roundUpTo((numBits * factor + 1) * sizeof(Count), CACHELINE_SIZE_BYTES);
    }

    [[nodiscard]] static constexpr Index requiredAlignment() noexcept { return alignof(Count); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept {
        if constexpr (Ops == Operations::RANK_ONLY) {
            return ranks.numT - 1;
        } else {
            return selectAnswers.numT;
        }
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR bool getBit(Index i) const noexcept {
        return rankOneUnchecked(i + 1) - this->rankOne(i);
    }

    // ** low-level functions **

    /// \brief Set the bit at position i. Assumes that this functions gets called for all bits in ascending order.
    /// \param i The ith bit will be set to newVal
    /// \param newVal new value of the ith bit.
    ADS_CPP20_CONSTEXPR void setBit(Index i, bool newVal = true) const noexcept {
        [[maybe_unused]] Index rankOfI = -1;
        if (i == 0) [[unlikely]] { // special case to allow iterating over and setting bits multiple times
            ranks.ptr[0] = 0;      // reset numZeros() to forget any previous calls to setBit()
        }
        if constexpr (Ops == Operations::SELECT_ONLY) { // ranks[0] counts the number of ones
            rankOfI = ranks.ptr[0];
            ranks.ptr[0] += newVal;
        } else {
            ranks.ptr[i + 1] = ranks[i] + newVal;
        }
        if constexpr (Ops == Operations::BOTH) {
            rankOfI = ranks[i];
        }
        if constexpr (Ops != Operations::RANK_ONLY) {
            if (newVal) {
                selectAnswers.ptr[rankOfI] = i;
            } else {
                Index rankZero = i - rankOfI;
                selectAnswers.ptr[size() - 1 - rankZero] = i;
            }
        }
    }

    ADS_CPP20_CONSTEXPR void buildMetadata() noexcept {
        if constexpr (Ops != Operations::SELECT_ONLY) {
            ranks.ptr[0] = 0; // may or may not have already been set
        }
    }

    // ** bitvector operations **

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numOnes() const noexcept {
        if constexpr (Ops == Operations::SELECT_ONLY) {
            return ranks[0];
        } else {
            return rankOneUnchecked(size());
        }
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOneUnchecked(Index i) const noexcept {
        if constexpr (Ops == Operations::SELECT_ONLY) {
            return this->template rankFallback<true>(i); // there is no fundamentally better solution
        } else {
            return ranks[i];
        }
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectOneUnchecked(Index i) const noexcept {
        if constexpr (Ops == Operations::RANK_ONLY) {
            return this->template inefficientSelect<true>(i); // there is no fundamentally better solution
        } else {
            return selectAnswers[i];
        }
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectZeroUnchecked(Index i) const noexcept {
        if constexpr (Ops == Operations::RANK_ONLY) {
            return this->template selectFallback<false>(i); // there is no fundamentally better solution
        } else {
            return selectAnswers[size() - 1 - i];
        }
    }
};


#ifdef ADS_HAS_CPP20
static_assert(IsBitvec<TrivialBitvec<>>);
#endif
} // namespace ads

#endif // ADS_TRIVIAL_BITVEC_HPP
