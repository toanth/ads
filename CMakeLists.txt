cmake_minimum_required(VERSION 3.25)
project(bitvector)

set(CMAKE_CXX_STANDARD 17)

add_executable(bitvector src/main.cpp include/bitvector.hpp include/common.hpp include/bitvec_layouts.hpp)

add_subdirectory(tests)
