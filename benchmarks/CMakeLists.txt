cmake_minimum_required(VERSION 3.25)

project(benchmarks)

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


#if (MSVC)
#    add_compile_options(/W4 /WX)
#else ()
#    add_compile_options(-Wall -Wextra -Wpedantic -Wno-sign-compare) # TODO: Enable Wsign-compare
#endif ()

if (MSVC)
    #    TODO
else ()
    add_compile_options(-fno-finite-math-only)
endif ()

add_executable(benchmarks bitvector_benchmarks.cpp ../include/bit.hpp eliasfano_benchmarks.cpp eliasfano_benchmarks.cpp benchmarks_common.hpp ../include/bitvector/bitvec_base.hpp ../include/bitvector/cache_efficient_rank_bitvec.hpp ../include/bitvector/efficient_rank_bitvec.hpp)

target_link_libraries(benchmarks benchmark::benchmark benchmark::benchmark_main)

target_include_directories(benchmarks PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../)
