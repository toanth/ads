#ifndef ADS_BITVECTOR_HPP
#define ADS_BITVECTOR_HPP

#include "bitvec_layouts.hpp"
#include "common.hpp"
#include <cassert>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <iostream>// TODO: Remove
#include <memory>
#include <stdexcept>

namespace ads {


template<ADS_LAYOUT_CONCEPT Layout>
class Bitvector : private Layout {
    using Base = Layout;
    using Base::getBlockCount;
    using Base::getElem;
    using Base::getSuperblockCount;
    using Base::setBlockCount;
    Index numBits;

    Bitvector(Index numBits, Index numElems) noexcept : Base(numElems), numBits(numBits) {
    }

public:
    using Base::superblockSize;

    explicit Bitvector(Index numBits) noexcept : Bitvector(numBits, ((numBits + 63) / 64 + superblockSize() - 1) / superblockSize() * superblockSize()) {}

    explicit Bitvector(std::string_view str, Index base = 2) : Bitvector(Index(str.size()) * log2(Elem(base))) {
        if (base != 2 && base != 4 && base != 16) [[unlikely]] {
            throw std::invalid_argument("base must be one of 2, 4, or 16");
        }
        const Index log2ofBase = log2(base);
        const std::size_t charsPerElem = 64 / log2ofBase;
        Index i = 0;
        for (Index superblock = 0; superblock < numSuperblocks(); ++superblock) {
            for (Index j = 0; j < superblockSize(); ++j) {
                if (str.empty()) break;
                std::size_t numToParse = std::min(charsPerElem, str.size());
                std::string_view toParse = str.substr(0, numToParse);
                std::uint64_t res;
                auto err = std::from_chars(toParse.data(), toParse.data() + toParse.size(), res, base);
                if (err.ec != std::errc()) [[unlikely]] {
                    throw std::invalid_argument(std::make_error_code(err.ec).message());
                } else if (err.ptr != toParse.data() + toParse.size()) [[unlikely]] {
                    throw std::invalid_argument("invalid character found");
                }
                if (numToParse < charsPerElem) {
                    res <<= (charsPerElem - numToParse) * log2ofBase;
                }
                str.remove_prefix(numToParse);
                getElem(i) = res;
                ++i;
            }
            buildMetadata(superblock);// TODO: Compute metadata while iterating over elements for cache efficiency
        }
        assert(str.empty());
    }


    void buildMetadata(Index superblockIdx) noexcept {
        assert(superblockIdx >= 0 && superblockIdx < numSuperblocks());
        if (superblockIdx == 0) {
            getSuperblockCount(0) = 0;
        }
        auto s = superblockElems(superblockIdx);
        Index inSuperblockSoFar = 0;
        for (Index i = 0; i < superblockSize(); ++i) {
            //            std::cout << "setting block count for block " << i << " to " << inSuperblockSoFar << ", s[i] is "
            //                      << std::hex << s[i] << std::dec << ", popcount " << popcount(s[i]) << std::endl;
            setBlockCount(superblockIdx, i, inSuperblockSoFar);
            inSuperblockSoFar += popcount(s[i]);
        }
        assert(inSuperblockSoFar <= superblockSize() * 64);
        if (superblockIdx < numSuperblocks() - 1) {
            getSuperblockCount(superblockIdx + 1) = getSuperblockCount(superblockIdx) + inSuperblockSoFar;
        }
    }

    [[nodiscard]] Index rankOne(Index pos) const {
        if (pos >= sizeInBits() || pos < 0) [[unlikely]] {
            throw std::invalid_argument("invalid position for rank query");
        }
        if (pos == 0) [[unlikely]] { return 0; }
        --pos;
        Index elemIdx = pos / 64;
        Elem mask = Elem(-1) << (63 - pos % 64);
        Index superblockIdx = elemIdx / superblockSize();
        Index blockIdx = elemIdx % superblockSize();
        return getSuperblockCount(superblockIdx) + getBlockCount(superblockIdx, blockIdx) + popcount(element(elemIdx) & mask);
    }

    [[nodiscard]] Index rankZero(Index pos) const {
        return pos - rankOne(pos);
    }

    //    explicit Bitvector(const char* name) noexcept {
    //        // use c io because that's faster than <iostreams>
    //        FILE* file = std::fopen(name, "r");
    //        if (!file) [[unlikely]]
    //            throw std::invalid_argument("couldn't open file");
    //        // TODO: posix_fadvise?
    //        constexpr static Index bufferSize = 1024 * 16 / sizeof(Elem);
    //        Elem buffer[bufferSize + 1];
    //        std::fread(buffer, sizeof(Elem), 1, file);
    //        Index numElems = Index(buffer[0]);// the number of entries in the bit vector
    //        assert(numElems >= 0 && std::size_t(numElems) == buffer[0]);
    //        vec = makeUniqueForOverwrite<Elem>(numElems);
    //        blocks = vec.get() + numElems;
    //        superBlocks = blocks + numBlocks(numElems);
    //        while (remaining > 0) {
    //            std::size_t numRead = std::fread(buffer, sizeof(Elem), std::min(bufferSize, remaining), file);
    //            assert(numRead < remaining);
    //            if (numRead < bufferSize && remaining != nunRead) [[unlikely]] {
    //                throw std::invalid_argument("io error");
    //            }
    //            remaining -= numRead;
    //        }
    //    }


    // Idea: superblock size is 4 Elems = 32 Byte = 256 bit, use 1 Byte to store number of 0s in superblock -- numbers in range [0, 8 * 24]
    // with s' = 256, we have s = 2^4 = 16 = (log n) / 2, therefore n = 2 ^ 32
    // n + 8 + n / 8 <= 64 <=> n <= 56 * 8 / 9 => 48, ergo 48 + 8 + 6 bytes, 48 Byte =
    // cache efficient layout: superblock bitvector (32 Byte) + Superblock count (8 Byte) + block count (1 Byte) * 4 == 44 Byte, wasting 20 Byte
    //or : Superblock size is less than 2^16 bit = 65636 bit = 4096 Byte = 512 Elems because s = (log n) / 2 <= 32, s' = s * s <= 1024

    [[nodiscard]] Index sizeInBits() const noexcept {
        return numBits;
    }

    [[nodiscard]] Index sizeInElems() const noexcept {
        return (numBits + 63) / 64;
    }

    [[nodiscard]] Index numSuperblocks() const noexcept {
        return (sizeInElems() + superblockSize() - 1) / superblockSize();
    }

    [[nodiscard]] Span<Elem> superblockElems(Index superblockIdx) noexcept {
        return Span<Elem>(&getElem(superblockIdx * superblockSize()), superblockSize());
    }
    [[nodiscard]] Span<const Elem> superblockElems(Index superblockIdx) const noexcept {
        return Span<const Elem>(&getElem(superblockIdx * superblockSize()), superblockSize());
    }


    [[nodiscard]] Elem element(Index elemIdx) const noexcept {
        return getElem(elemIdx);
    }

    [[nodiscard]] bool bit(Index bitIndex) const noexcept {
        return (element(bitIndex / 64) >> (63 - bitIndex % 64)) & 0x1;
    }
};

}// namespace ads

#endif// ADS_BITVECTOR_HPP
