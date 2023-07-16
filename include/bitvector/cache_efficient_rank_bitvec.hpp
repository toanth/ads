#ifndef ADS_CACHE_EFFICIENT_RANK_BITVEC_HPP
#define ADS_CACHE_EFFICIENT_RANK_BITVEC_HPP

#include "base/superblock_bitvec.hpp"

namespace ads {

// This layout wastes a lot of memory; it needs roughly double the total memory of the actual bit vector
// Still, for small bitvectors, this may be a good option
class [[nodiscard]] CacheEfficientRankBitvec : public SuperblockBitvecBase<CacheEfficientRankBitvec, 4, 1, Limb, 4> {

    using Base = SuperblockBitvecBase<CacheEfficientRankBitvec, 4, 1, Limb, 4>;
    friend Base;
    friend Base::Base;
    friend Base::Base::Base; // :D

    Array<CacheLine> vec = Array<CacheLine>();

    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<const CacheLine> getCacheLineArray() const noexcept {
        return {vec.ptr, vec.numT};
    }
    [[nodiscard]] ADS_CPP20_CONSTEXPR Span<CacheLine> getCacheLineArray() noexcept { return {vec.ptr, vec.numT}; }

    explicit ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(UninitializedTag, Index numBits, CacheLine* mem = nullptr) noexcept
        : Base(numBits, mem), vec(this->allocation.memory(), this->allocation.sizeInTs()) {
        [[maybe_unused]] auto ptr = allocation.memory();
        ADS_ASSUME_ALIGNED(ptr, CACHELINE_SIZE_BYTES);
        ADS_ASSUME(this->allocation.isEnd(vec.end()));
    }

public:
    using typename Base::BlockRank;
    using typename Base::SuperblockRank;

    constexpr CacheEfficientRankBitvec() noexcept = default;

    ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(Index numBits, Limb fill, CacheLine* mem = nullptr)
        : CacheEfficientRankBitvec(UninitializedTag{}, numBits, mem) {
        this->fill(fill);
    }

    explicit ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(Span<const Limb> limbs, CacheLine* mem = nullptr) noexcept
        : CacheEfficientRankBitvec(limbs, limbs.size() * 64, mem) {}

    ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(Span<const Limb> limbs, Index numBits, CacheLine* mem = nullptr) noexcept
        : CacheEfficientRankBitvec(UninitializedTag{}, numBits, mem) {
        this->copyFrom(limbs);
    }

    explicit ADS_CPP20_CONSTEXPR CacheEfficientRankBitvec(std::string_view str, Index base = 2, CacheLine* mem = nullptr) noexcept
        : CacheEfficientRankBitvec(UninitializedTag{}, str.size() * intLog2(base), mem) {
        initFromStr(str, base);
    }

    //    [[nodiscard]] constexpr static Index numLimbsInSuperblock() noexcept { return 4; }

    //    [[nodiscard]] constexpr static Index numLimbsInBlock() noexcept { return 1; }

    [[nodiscard]] static constexpr Index allocatedSizeInBytesForBits(Index numBits) noexcept {
        return (1 + roundUpDiv(numBits, superblockSize())) * 64;
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR SuperblockRank getSuperblockRank(Index i) const noexcept {
        ADS_ASSUME(i >= 0);
        // the number of superblock ranks is one larger than the number of superblocks
        ADS_ASSUME(i <= derived().numSuperblocks());
        return vec[i][numLimbsInSuperblock()];
    }

    ADS_CPP20_CONSTEXPR void setSuperblockRank(Index i, SuperblockRank newVal) noexcept {
        ADS_ASSUME(i >= 0);
        // the number of superblock ranks is one larger than the number of superblocks
        ADS_ASSUME(i <= derived().numSuperblocks());
        vec.ptr[i][numLimbsInSuperblock()] = newVal;
    }


    ADS_CPP20_CONSTEXPR void setBlockRank(Index i, BlockRank newVal) noexcept {
        ADS_ASSUME(i >= 0);
        ADS_ASSUME(i < derived().numBlocks());
        setBlockRank_(i / numBlocksInSuperblock(), i % numBlocksInSuperblock(), newVal);
    }

    ADS_CPP20_CONSTEXPR void setBlockRank_(Index superblockIdx, Index blockInSuperblock, BlockRank newVal) noexcept {
        ADS_ASSUME(superblockIdx >= 0);
        ADS_ASSUME(superblockIdx < derived().numSuperblocks());
        ADS_ASSUME(blockInSuperblock >= 0);
        ADS_ASSUME(blockInSuperblock < derived().numBlocksInSuperblock());
        BitwiseAccess<8>::setBits(getCompleteCacheLine(superblockIdx).limbs, numLimbsInSuperblock() + 1, blockInSuperblock, newVal);
    }

    [[nodiscard]] ADS_CPP20_CONSTEXPR Index getBlockRank(Index blockIdx) const noexcept {
        ADS_ASSUME(blockIdx >= 0);
        ADS_ASSUME(blockIdx < numAccessibleBlocks());
        return Index(BitwiseAccess<8>::getBits(getCompleteCacheLine(blockIdx / numBlocksInSuperblock()).limbs,
                numLimbsInSuperblock() + 1, blockIdx % numBlocksInSuperblock()));
    }
};

#ifdef ADS_HAS_CPP20
static_assert(IsNormalBitvec<CacheEfficientRankBitvec>);
#endif

} // namespace ads

#endif // ADS_CACHE_EFFICIENT_RANK_BITVEC_HPP
