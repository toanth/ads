cmake_minimum_required(VERSION 3.25)

project(project)


if (MSVC)
    #    TODO
else ()
    add_compile_options(-march=native -flto=auto)
    add_link_options(-flto=auto)
endif ()

add_compile_options(project -DADS_NO_FALLBACKS)

add_executable(project main.cpp ../include/bitvector/bitvec_base.hpp ../include/bitvector/cache_efficient_rank_bitvec.hpp ../include/bitvector/efficient_rank_bitvec.hpp ../include/bitvector/trivial_bitvec.hpp ../include/bitvector/recursive_bitvec.hpp)
