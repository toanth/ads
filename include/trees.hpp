//
// Created by tobias on 01/06/23.
//

#ifndef BITVECTOR_TREES_HPP
#define BITVECTOR_TREES_HPP

#include "bitvector.hpp"
namespace ads {

template<typename T>
Bitvector<> createDfudsOf2dMinHeap(Span<const T> values) {
    Bitvector<> bv(2 * values.size());
    std::vector<Index> stack;
    Elem current = 0;
    for (Index i = values.size() - 1, elemIdx = (values.size() - 1) / 64, bitIdx = 63 - (values.size() - 1) % 64; i + 1 > 0; --i, ++bitIdx) {
        if (bitIdx == 64) {
            bv.setElem(elemIdx, current);
            bitIdx = 0;
            --elemIdx;
            current = 0;
        }
        current |= (Elem(1) << bitIdx);
    }
}
}// namespace ads


#endif//BITVECTOR_TREES_HPP
