cmake_minimum_required(VERSION 3.25)
project(bitvector)

set(CMAKE_CXX_STANDARD 17)

add_executable(bitvector main.cpp bitvector.hpp common.hpp bitvec_layouts.hpp)
