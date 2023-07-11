#ifndef ADS_CONCEPT_HPP
#define ADS_CONCEPT_HPP

#include "../../bit.hpp"
#include "../../bit_access.hpp"
#include "../../common.hpp"
#include "../../rand_access_iter.hpp"

namespace ads {

/// Note that the following concepts also apply in C++17 mode, they are just not checked by the implementation.
#ifdef ADS_HAS_CPP20
/// Bitvectors organize their bit sequence in groups of different sizes, where each group contains one or more of the
/// preceding groups, and each instance of a group has the same size (with the possible exception of the last
/// superblock). The groups are: Bit, Byte, Limb, Block, Cache Line, Superblock. In general, bitvectors only have the
/// "Bit" group, but a `NormalBitvec` also consists of Bytes, Blocks and Cache Lines, while a SuperblockBitvec also has
/// superblocks. A Bit is (obviously) a group of 1 bit, a Byte a group of 8 bit, a Limb a group of 64 bit, a block a
/// group of implementation-defined size that is a multiple of 64 and a divisor of bv.numBitsStoredPerCacheLine, which
/// is the number of bits the bitvector bv stores in a Cache Line, while a superblock is an implementation-defined
/// multiple of that number. Note that the number of bits stored in a cache line may be less than CACHELINE_SIZE_BYTES *
/// 8 because a bitvector may store meta data in the same cache line as the bit sequence. Blocks and CacheLines must be
/// a sequence of Limbs contiguous in memory, a Cache Line must be contained within and aligned to an actual cache line;
/// a Limb must be represented by the type Limb (for non-NormalBitvecs, all of this is irrelevant as they don't have
/// those groups). Because the number of bits is not necessarily a multiple of these groups sizes, bitvectors can extend
/// the stored bit sequence and provide the numAccessible<Group> functions in addition to the sizeIn<Group> functions.
enum class Group { Bit, Byte, Limb, Block, CacheLine, Superblock };

template<Group G>
using TypeFor = std::conditional_t<G == Group::Byte, Byte,
        std::conditional_t<G == Group::Limb, Limb, std::conditional_t<G == Group::CacheLine, CacheLine, void>>>;


/// This concept describes a very general bitvector, all implementations except for TrivialBitvec additionally model
/// IsNormalBitvec. Inheriting from BitvecBase provides default implementations for most required methods.
template<typename T>
concept IsBitvec = requires(T& t, const T& ct) {
    T();
    T(Index(), Limb());
    T(Index(), Limb(), (CacheLine*)nullptr);
    T("");
    T("", Index());
    T("", Index(), (CacheLine*)nullptr);

    // These static methods are called to reserve the necessary memory beforehand.
    // The std::integral_constant enforces that this method is constexpr.
    { T::template allocatedSizeForBitsIn<Group::Byte>(Index()) } -> std::convertible_to<Index>;
    { T::template allocatedSizeForBitsIn<Group::Limb>(Index()) } -> std::convertible_to<Index>;
    { T::template allocatedSizeForBitsIn<Group::CacheLine>(Index()) } -> std::convertible_to<Index>;

    std::integral_constant<Index, (T::requiredAlignment())>{};

    { ct.template allocatedSizeIn<Group::Byte>() } -> std::convertible_to<Index>;
    { ct.template allocatedSizeIn<Group::Limb>() } -> std::convertible_to<Index>;
    { ct.template allocatedSizeIn<Group::CacheLine>() } -> std::convertible_to<Index>;

    std::integral_constant<Index, T::requiredAlignment()>{};
    { T::allocatedSizeInBytesForBits(Index()) } -> std::convertible_to<Index>;
    { ct.allocatedSizeInBytes() } -> std::convertible_to<Index>;
    //    { ct.allocatedSizeInCacheLines() } -> std::convertible_to<Index>;

    { ct.size() } -> std::convertible_to<Index>;
    { ct.numBits() } -> std::convertible_to<Index>;

    { ct.template get<Group::Bit>(Index()) } -> std::convertible_to<bool>;

    { ct.getBit(Index()) } -> std::convertible_to<bool>;
    t.buildMetadata();


    { ct.numOnes() } -> std::convertible_to<Index>;
    { ct.numZeros() } -> std::convertible_to<Index>;
    { ct.template numBitsEqualTo<true>() } -> std::convertible_to<Index>;
    { ct.template numBitsEqualTo<false>() } -> std::convertible_to<Index>;

    { ct.rankOne(Index()) } -> std::convertible_to<Index>;
    { ct.rankZero(Index()) } -> std::convertible_to<Index>;
    { ct.template rank<true>(Index()) } -> std::convertible_to<Index>;
    { ct.template rank<false>(Index()) } -> std::convertible_to<Index>;
    { ct.selectOne(Index()) } -> std::convertible_to<Index>;
    { ct.selectZero(Index()) } -> std::convertible_to<Index>;
    { ct.template select<true>(Index()) } -> std::convertible_to<Index>;
    { ct.template select<false>(Index()) } -> std::convertible_to<Index>;
};


/// This concepts describes a bitvector that stores a sequence of bits, organized as an array of cache lines (which are
/// themselves arrays of limbs), where each cache line is split into a fixed number of fixed-sized blocks. The bitvector
/// may also store additional metadata in an unspecified format. Although all the following methods must be present,
/// inheriting from NormalBitvecBase provides default implementations for most.
template<typename T>
concept IsNormalBitvec = IsBitvec<T> && requires(T& t, const T& ct) {
    std::integral_constant<Index, T::blockSize()>{};
    std::integral_constant<Index, T::numLimbsInBlock()>{};
    std::integral_constant<Index, T::numBlocksInCacheLine()>{};
    std::integral_constant<Index, T::numLimbsInCacheLine()>{};
    std::integral_constant<Index, T::bytesPerBlockRank()>{};

    //    requires T::template numBitsIn<Group::Bit>() == 1;
    //    requires T::template numBitsIn<Group::Byte>() == 8;
    //    requires T::template numBitsIn<Group::Limb>() == 64;
    //    requires T::template numBitsIn<Group::Block>() % 64 == 0;
    //    requires T::template numBitsIn<Group::CacheLine>() <= CACHELINE_SIZE_BYTES * 8;
    //    requires T::template numBitsIn<Group::CacheLine>() % T::template numBitsIn<Group::Block>() == 0;

    ct.template num<Group::Bit>();
    ct.template num<Group::Byte>();
    ct.template num<Group::Limb>();
    ct.template num<Group::Block>();
    ct.template num<Group::CacheLine>();


    ct.template numAccessible<Group::Bit>();
    ct.template numAccessible<Group::Byte>();
    ct.template numAccessible<Group::Limb>();
    ct.template numAccessible<Group::Block>();
    ct.template numAccessible<Group::CacheLine>();

    std::integral_constant<Index, T::numAccessibleBytesForBits(Index())>{};

    { ct.numLimbs() } -> std::convertible_to<Index>;
    { ct.numAccessibleLimbs() } -> std::convertible_to<Index>;
    { ct.numAccessibleCacheLines() } -> std::convertible_to<Index>;
    { ct.numAccessibleBlocks() } -> std::convertible_to<Index>;


    { ct.template get<Group::Byte>(Index()) } -> std::convertible_to<Byte>;
    { ct.template get<Group::Limb>(Index()) } -> std::convertible_to<Limb>;
    ct.template get<Group::Block>(Index());
    t.template get<Group::Block>(Index());
    ct.template get<Group::CacheLine>(Index());
    t.template get<Group::CacheLine>(Index());

    { ct.numBlocks() } -> std::convertible_to<Index>;
    { ct.getLimb(Index()) } -> std::convertible_to<const U64&>;
    { ct.numCacheLines() } -> std::convertible_to<Index>;
    ct.getCacheLine(Index());
    //    { ct.numLimbsStoredPerCacheLine() } -> std::convertible_to<Index>; // TODO: Remove

    t.setLimb(Index(), Limb());
    t.setBit(Index());
    t.setBit(Index(), bool());
};

/// This concepts describes a bitvector that stores a sequence of bits, organized as an array of limbs, and additionally
/// partitions them in blocks and superblocks (where a superblock may consists of only 1 block) of fixed size to answer
/// rank queries.
/// Although all the following methods must be present, inheriting from SuperblockBitvecBase provides default implementations for most.
template<typename T>
concept IsSuperblockBitvec = IsNormalBitvec<T> && requires(T& t, const T& ct) {
    typename T::BlockRank;
    typename T::SuperblockRank;

    //    requires T::template numBitsIn<Group::Superblock>() % T::template numBitsIn<Group::Block>() == 0;

    std::integral_constant<Index, T::superblockSize()>{};
    std::integral_constant<Index, T::numLimbsInSuperblock()>{};
    std::integral_constant<Index, T::numCacheLinesInSuperblock()>{};
    std::integral_constant<Index, T::bytesPerSuperblockRank()>{};
    std::integral_constant<Index, T::numBlocksInSuperblock()>{};
    std::integral_constant<Index, T::numSuperblocksForBits(Index())>{};

    { ct.getSuperblockRank(Index()) } -> std::convertible_to<const U64&>;
    t.setSuperblockRank(Index(), Index());
    { ct.getBlockRank(Index()) } -> std::convertible_to<Index>;
    t.setBlockRank(Index(), Index());
    t.setBlockRank_(Index(), Index(), Index());
    { ct.numSuperblocks() } -> std::convertible_to<Index>;
};


#define ADS_BITVEC_CONCEPT IsBitvec
#define ADS_NORMAL_BITVEC_CONCEPT IsNormalBitvec
#define ADS_SUPERBLOCK_BITVEC_CONCEPT IsSuperblockBitvec

#else // ADS_HAS_CPP20

#define ADS_BITVEC_CONCEPT class
#define ADS_NORMAL_BITVEC_CONCEPT class
#define ADS_SUPERBLOCK_BITVEC_CONCEPT class

namespace detail {

template<typename T, typename = void>
struct IsBitvecImpl : std::false_type {};

template<typename T>
struct IsBitvecImpl<T, std::void_t<decltype(T::getBit(1))>> : std::true_type {};

template<typename T, typename = void>
struct IsNormalBitvecImpl : std::false_type {};

template<typename T>
struct IsNormalBitvecImpl<T, std::void_t<decltype(T().getLimb(1))>> : std::true_type {};

template<typename T, typename = void>
struct IsSuperblockBitvecImpl : std::false_type {};

template<typename T>
struct IsSuperblockBitvecImpl<T, std::void_t<decltype(T::numSuperblocksForBits(1))>> : std::true_type {};

}; // namespace detail

template<typename T>
constexpr static bool IsBitvec = detail::IsNormalBitvecImpl<T>::value;

template<typename Bitvec>
constexpr static bool IsNormalBitvec = IsBitvec<Bitvec> && detail::IsNormalBitvecImpl<T>::value;

template<typename Bitvec>
constexpr static bool IsSuperblockBitvec = IsNormalBitvec<Bitvec> && detail::IsSuperblockBitvecImpl<T>::value;

#endif // ADS_HAS_CPP20


/// These lambdas are helpful to implement views over the bitvectors
template<Group G>
[[maybe_unused]] constexpr static inline auto getGroupFunc
        = [](const auto& bv, Index i) -> bool { return bv.template get<G>(i); };

template<bool IsOne>
[[maybe_unused]] constexpr static inline auto rankFunc
        = [](const auto& bv, Index i) -> Index { return bv.template rank<IsOne>(i); };

template<bool IsOne>
[[maybe_unused]] constexpr static inline auto selectFunc
        = [](const auto& bv, Index i) -> Index { return bv.template select<IsOne>(i); };


/// Comparison operators

template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] ADS_CPP20_CONSTEXPR bool operator==(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return lhs.numBits() == rhs.numBits() && std::equal(lhs.bitView().begin(), lhs.bitView().end(), rhs.bitView().begin());
}
template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] ADS_CPP20_CONSTEXPR bool operator!=(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return !(rhs == lhs);
}

#ifdef ADS_HAS_CPP20

template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] constexpr std::strong_ordering operator<=>(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return lexicographicalCompareThreeWay(
            lhs.bitView().begin(), lhs.bitView().end(), rhs.bitView().begin(), rhs.bitView().end());
}

#else

template<ADS_BITVEC_CONCEPT Bitvec1, ADS_BITVEC_CONCEPT Bitvec2>
[[nodiscard]] bool operator<(const Bitvec1& lhs, const Bitvec2& rhs) noexcept {
    return std::lexicographical_compare(lhs.bitView().begin(), lhs.bitView().end(), rhs.bitView().begin(), rhs.bitView().end());
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

#endif // ADS_CONCEPT_HPP
