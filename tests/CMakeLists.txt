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
    add_compile_options(-Wall -Wextra -Wpedantic -g)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-fconstexpr-steps=10000000)
    endif ()
endif ()

add_executable(tests common_tests.cpp bit_tests.cpp bitvector_tests.cpp eliasfano_tests.cpp rmq_tests.cpp ../include/succinct_rmq.hpp rmmtree_tests.cpp ../include/bitvector/select_bitvec.hpp)

target_link_libraries(tests gtest gtest_main)
