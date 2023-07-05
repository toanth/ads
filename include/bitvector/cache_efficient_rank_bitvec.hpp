#ifndef ADS_CACHE_EFFICIENT_RANK_BITVEC_HPP
#define ADS_CACHE_EFFICIENT_RANK_BITVEC_HPP

#include "bitvec_base.hpp"

namespace ads {

// This layout wastes a lot of memory; it needs roughly double the total memory of the actual bit vector
// Still, for small bitvectors, this may be a good option
class [[nodiscard]] CacheEfficientRankBitvec : public RankBitvecBase<CacheEfficientRankBitvec, std::uint8_t> {

    using Base = RankBitvecBase<CacheEfficientRankBitvec, std::uint8_t>;
    friend Base;
    friend Base::Base;
    friend Base::Base::Base; // :D

    View<Limb> vec = View<Limb>();
    Index numBits = 0;

    [[nodiscard]] static constexpr Index limbIdxInArray(Index i) noexcept {
        ADS_ASSUME(i >= 0);
        return (i / numLimbsInSuperblock()) * (U64_PER_CACHELINE) + i % numLimbsInSuperblock();
    }

    [[nodiscard]] static constexpr Index superblockRankIdxInArray(Index superblockIdx) noexcept {
        ADS_ASSUME(superblockIdx >= 0);
        return superblockIdx * U64_PER_CACHELINE + numLimbsInSuperblock();
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const Limb* getLimbArray() const noexcept { return vec.ptr; }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Limb* getLimbArray() noexcept { return vec.ptr; }

    [[nodiscard]] ADS_CPP20_CONSTEXPR const View<SuperblockRank>& getSuperblockRankArray() const noexcept {
        return vec;
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR View<SuperblockRank>& getSuperblockRankArray() noexcept { return vec; }

    explicit ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(UninitializedTag, Index numBits, Limb* mem = nullptr) noexcept
        : Base(numBits, mem), vec(this->allocation.memory(), this->allocation.size()), numBits(numBits) {}

public:
    using typename Base::BlockRank;
    using typename Base::SuperblockRank;

    constexpr CacheEfficientRankBitvec() noexcept = default;

    ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(Index numBits, Elem fill, Elem* mem = nullptr)
        : CacheEfficientRankBitvec(UninitializedTag{}, numBits, mem) {
        this->fill(fill);
    }

    explicit ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(Span<const Limb> limbs, Limb* mem = nullptr) noexcept
        : CacheEfficientRankBitvec(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(Span<const Limb> limbs, Index numBits, Limb* mem = nullptr) noexcept
        : CacheEfficientRankBitvec(UninitializedTag{}, numBits, mem) {
        this->copyFrom(limbs);
    }

    explicit ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(std::string_view str, Index base = 2, Limb* mem = nullptr) noexcept
        : CacheEfficientRankBitvec(UninitializedTag{}, str.size() * intLog2(base), mem) {
        initFromStr(str, base);
    }

    [[nodiscard]] constexpr static Index numLimbsInSuperblock() noexcept { return 4; }

    [[nodiscard]] constexpr static Index numLimbsInBlock() noexcept { return 1; }

    [[nodiscard]] static constexpr Index allocatedSizeInLimbsForBits(Index numBits) noexcept {
        return (1 + roundUpDiv(numBits, superblockSize())) * U64_PER_CACHELINE;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index size() const noexcept { return numBits; }


    [[nodiscard]] ADS_CPP20_CONSTEXPR Index getBlockRank(Index blockIdx) const noexcept {
        Index superblockIdx = blockIdx / numBlocksInSuperblock();
        blockIdx %= numBlocksInSuperblock();
        ADS_ASSUME(superblockIdx >= 0);
        ADS_ASSUME(superblockIdx < numSuperblocks());
        return (Index)BitwiseAccess<8>::getBits(vec.ptr, superblockRankIdxInArray(superblockIdx) + 1, blockIdx);
    }

    ADS_CPP20_CONSTEXPR void setBlockRank_(Index superblockIdx, Index blockIdx, Index value) noexcept {
        ADS_ASSUME(superblockIdx >= 0 && superblockIdx < numSuperblocks() && blockIdx >= 0 && blockIdx < numBlocks());
        ADS_ASSUME(value >= 0);
        Index i = superblockRankIdxInArray(superblockIdx) + 1;
        BitwiseAccess<8>::setBits(vec.ptr, i, blockIdx, value);
    }
};

#ifdef ADS_HAS_CPP20
static_assert(IsNormalBitvec<CacheEfficientRankBitvec>);
#endif

} // namespace ads

#endif // ADS_CACHE_EFFICIENT_RANK_BITVEC_HPP
