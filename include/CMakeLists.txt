cmake_minimum_required(VERSION 3.25)
project(bitvector_lib)

set(CMAKE_CXX_STANDARD 17)

add_library(bitvector_lib include/bitvector.hpp include/common.hpp include/bitvec_layouts.hpp)
