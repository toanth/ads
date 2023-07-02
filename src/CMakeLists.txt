cmake_minimum_required(VERSION 3.25)

project(project)


if (MSVC)
    #    TODO
else ()
    add_compile_options(-Ofast -march=native)
endif ()

add_executable(project main.cpp ../include/bitvector/bitvec_base.hpp ../include/bitvector/cache_efficient_rank_bitvec.hpp ../include/bitvector/efficient_rank_bitvec.hpp)
