cmake_minimum_required(VERSION 3.25)
project(bitvector)

set(CMAKE_CXX_STANDARD 17)

add_executable(bitvector src/main.cpp include/bitvector.hpp include/common.hpp include/bitvec_layouts.hpp include/elias_fano.hpp include/rand_access_iter.hpp include/naive_rmq.hpp include/linear_space_rmq.hpp include/succinct_rmq.hpp include/trees.hpp include/bit_access.hpp)

add_subdirectory(tests)
