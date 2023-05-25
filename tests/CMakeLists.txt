project(tests)
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(tests bitvector_tests.cpp ../include/elias_fano.hpp)

#target_link_libraries(tests bitvector)

target_link_libraries(tests gtest gtest_main)
