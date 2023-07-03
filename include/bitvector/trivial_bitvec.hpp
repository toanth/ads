#ifndef ADS_TRIVIAL_BITVEC_HPP
#define ADS_TRIVIAL_BITVEC_HPP

#include "../common.hpp"
#include "bitvec_base.hpp"

namespace ads {

/// \brief This class precomputes the answers to all possible rank and select queries and stores them in one large array.
/// `getBit(i)` is implemented as `rankOne(i+1) - rankOne(i)` for all i (there are size() + 1 rank entries).
/// For 32 bit rank and select indices, this uses 64 times as much space as the actual bitvector (which isn't stored),
/// so this class is only useful for small Bitvectors, such as in the second recursion level of the recursive bitvector
// TODO: Implement recursive bitvector
template<typename Count = std::uint32_t>
class TrivialBitvec : public BitvecBase<TrivialBitvec<Count>> {

    using Base = BitvecBase<TrivialBitvec<Count>>;
    friend Base;

    View<Count> ranks = View<Count>();
    View<Count> selectAnswers = View<Count>();

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
                ranks.ptr[bitNum] = previousRank;
                previousRank += bit;
                if (bit) {
                    selectAnswers.ptr[numOnesSoFar++] = bitNum;
                } else {
                    selectAnswers.ptr[size() - ++numZerosSoFar] = bitNum;
                }
            }
        }
        ranks.ptr[size()] = previousRank;
    }

    ADS_CPP20_CONSTEXPR TrivialBitvec(UninitializedTag, Index numBits, Limb* mem) noexcept
        : Base(numBits, mem), ranks(this->allocation.memory(), numBits + 1), selectAnswers(ranks.ptr + numBits + 1, numBits) {}

public:
    constexpr TrivialBitvec() noexcept = default;
    ADS_CPP20_CONSTEXPR TrivialBitvec(Index size, Limb fill, Limb* mem = nullptr) noexcept
        : TrivialBitvec(UninitializedTag{}, size, mem) {
        constructFromLimbRange(repeatView(fill, roundUpDiv(size, 64)));
    }

    explicit ADS_CPP20_CONSTEXPR TrivialBitvec(Span<const Limb> limbs, Limb* mem = nullptr) noexcept
        : TrivialBitvec(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR TrivialBitvec(Span<const Limb> limbs, Index numBits, Limb* mem = nullptr) noexcept
        : TrivialBitvec(UninitializedTag{}, numBits, mem) {
        constructFromLimbRange(limbs);
    }


    ADS_CPP20_CONSTEXPR TrivialBitvec(std::string_view str, Index base = 2, Limb* mem = nullptr)
        : TrivialBitvec(UninitializedTag{}, str.size() * intLog2(base), mem) {
        constructFromLimbRange(this->limbViewFromStringView(str, base));
    }

    [[nodiscard]] static ADS_CPP20_CONSTEXPR Index allocatedSizeInLimbsForBits(Index numBits) noexcept {
        return numBits * sizeof(Count) * 2 + 1;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return selectAnswers.numT; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR bool getBit(Index i) const noexcept { return rankOne(i + 1) - rankOne(i); }

    ADS_CPP20_CONSTEXPR void buildMetadata() noexcept {
        // do nothing
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index numOnes() const noexcept { return rankOne(size()); }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index rankOne(Index i) const noexcept { return ranks[i]; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectOne(Index i) const noexcept { return selectAnswers[i]; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index selectZero(Index i) const noexcept { return selectAnswers[size() - 1 - i]; }
};


#ifdef ADS_HAS_CPP20
static_assert(IsBitvec<TrivialBitvec<>>);
#endif
} // namespace ads

#endif // ADS_TRIVIAL_BITVEC_HPP
