
#ifndef BITVECTOR_SUCCINCT_RMQ_HPP
#define BITVECTOR_SUCCINCT_RMQ_HPP

#include "bitvector.hpp"
#include "common.hpp"
namespace ads {

class SuccinctRMQ {

    Bitvector<> bv;


public:
    constexpr static const char name[] = "Succinct Rmq";

    SuccinctRMQ() = default;
    template<typename T>
    explicit SuccinctRMQ(Span<const T>()) {
    }
    template<typename T>
    SuccinctRMQ(std::unique_ptr<T> ptr, Index length) : SuccinctRMQ(Span<const T>(ptr.get(), length)) {}
    template<typename T>
    SuccinctRMQ(std::initializer_list<T> list) : SuccinctRMQ(Span<const T>(list.begin(), list.end())) {}

    [[nodiscard]] Index rmq(Index lower, Index upper) const noexcept {
    }
    [[nodiscard]] Index operator()(Index lower, Index upper) const noexcept {
        return rmq(lower, upper);
    }

    [[nodiscard]] Index sizeInBits() const noexcept {}
};
};

}// namespace ads

#endif//BITVECTOR_SUCCINCT_RMQ_HPP
