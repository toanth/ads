#ifndef BITVECTOR_RAND_ACCESS_ITER_HPP
#define BITVECTOR_RAND_ACCESS_ITER_HPP

#include "common.hpp"

namespace ads {

template<typename Container, typename Projection, typename Ref = decltype(std::declval<Projection&>()(std::declval<Container&>(), Index()))>
class RandAccessIter {

    const Container* cPtr;
    [[no_unique_address]] Projection p;
    Index i;

public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::remove_reference_t<Ref>;// don't remove const
    using difference_type = Index;
    using pointer = value_type*;
    using reference = Ref;

    RandAccessIter() : cPtr(nullptr), p(Projection()), i(0) {}
    explicit RandAccessIter(const Container& c, Projection p, Index pos = 0) : cPtr(&c), p(p), i(pos) {}

    RandAccessIter(const RandAccessIter&) = default;
    RandAccessIter& operator=(const RandAccessIter& other) {
        cPtr = other.cPtr;
        i = other.i;
        return *this;
    }

    RandAccessIter& operator++() {
        return *this += 1;
    }
    RandAccessIter operator++(int) {
        RandAccessIter copy(*this);
        ++*this;
        return copy;
    }
    RandAccessIter& operator--() {
        return *this -= 1;
    }
    RandAccessIter operator--(int) {
        RandAccessIter copy(*this);
        --*this;
        return copy;
    }
    RandAccessIter& operator+=(Index n) {
        i += n;
        return *this;
    }
    RandAccessIter& operator-=(Index n) {
        i -= n;
        return *this;
    }

    friend RandAccessIter operator+(RandAccessIter iter, Index n) {
        return iter += n;
    }
    friend RandAccessIter operator+(Index n, RandAccessIter iter) {
        return iter += n;
    }
    friend RandAccessIter operator-(RandAccessIter iter, Index n) {
        return iter -= n;
    }
    friend Index operator-(RandAccessIter a, RandAccessIter b) {
        return a.i - b.i;
    }

    reference operator*() const {
        return operator[](0);
    }

    reference operator[](Index n) const {
        return p(*cPtr, i + n);
    }

    friend bool operator==(const RandAccessIter& lhs, const RandAccessIter& rhs) {
        return lhs.i == rhs.i;
    }
    friend bool operator!=(const RandAccessIter& lhs, const RandAccessIter& rhs) {
        return !(rhs == lhs);
    }

#ifdef ADS_HAS_CPP20
    friend std::strong_ordering operator<=>(RandAccessIter, RandAccessIter) = default;
#else
    friend bool operator<(RandAccessIter lhs, RandAccessIter rhs) {
        return lhs.i < rhs.i;
    }
    friend bool operator>(RandAccessIter lhs, RandAccessIter rhs) {
        return rhs < lhs;
    }
    friend bool operator<=(RandAccessIter lhs, RandAccessIter rhs) {
        return !(rhs < lhs);
    }
    friend bool operator>=(RandAccessIter lhs, RandAccessIter rhs) {
        return !(lhs < rhs);
    }
#endif
};

}// namespace ads

#endif//BITVECTOR_RAND_ACCESS_ITER_HPP
