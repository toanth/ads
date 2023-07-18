cmake_minimum_required(VERSION 3.25)

project(benchmarks)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(GOOGLETEST_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/googletest-src")

include(FetchContent)
set(FETCHCONTENT_QUIET off)
get_filename_component(fc_base "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies" REALPATH BASE_DIR "${CMAKE_BINARY_DIR}")
set(FETCHCONTENT_BASE_DIR ${fc_base})

FetchContent_Declare(
        benchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.8.0
)

FetchContent_GetProperties(benchmark)
if (NOT benchmark_POPULATED)
    FetchContent_Populate(benchmark)
    add_subdirectory(${benchmark_SOURCE_DIR} benchmark)
endif ()

if (MSVC)
    #    TODO
else ()
    add_compile_options(-fno-finite-math-only -flto=auto)
    add_link_options(-flto=auto)
endif ()

add_compile_options(benchmarks -DADS_NO_FALLBACKS)

add_executable(benchmarks bitvector_benchmarks.cpp eliasfano_benchmarks.cpp ../include/bitvector/classical_bitvec.hpp ../include/bitvector/select_bitvec.hpp rmq_benchmarks.cpp)

target_link_libraries(benchmarks benchmark::benchmark benchmark::benchmark_main)

#target_include_directories(benchmarks PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../)


add_executable(profiling profiling.cpp ../include/bitvector/base/concept.hpp ../include/bitvector/base/bitvec.hpp ../include/bitvector/base/normal_bitvec.hpp ../include/bitvector/base/superblock_bitvec.hpp)
target_link_libraries(profiling benchmark::benchmark)
#target_include_directories(profiling PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../)
