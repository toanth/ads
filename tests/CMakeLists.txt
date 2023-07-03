cmake_minimum_required(VERSION 3.25)

project(tests)

include(FetchContent)
set(FETCHCONTENT_QUIET off)
get_filename_component(fc_base "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies" REALPATH BASE_DIR "${CMAKE_BINARY_DIR}")
set(FETCHCONTENT_BASE_DIR ${fc_base})

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.13.0
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_GetProperties(googletest)
if (NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} googletest)
endif ()
enable_testing()


if (MSVC)
    add_compile_options(/W4 /WX)
else ()
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-sign-compare) # TODO: Enable Wsign-compare
endif ()

add_executable(tests common_tests.cpp bit_tests.cpp bitvector_tests.cpp ../include/elias_fano.hpp eliasfano_tests.cpp rmq_tests.cpp ../include/nlogn_rmq.hpp ../include/linear_space_rmq.hpp ../include/succinct_rmq.hpp rmmtree_tests.cpp ../include/bit.hpp ../include/bitvector/bitvec_base.hpp ../include/bitvector/cache_efficient_rank_bitvec.hpp ../include/bitvector/trivial_bitvec.hpp)

target_link_libraries(tests gtest gtest_main)
