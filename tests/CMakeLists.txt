project(tests)
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

if (MSVC)
    add_compile_options(/W4 /WX)
else ()
    #    add_compile_options(-Wall -Wextra -Wpedantic -Wno-sign-compare) # TODO: Enable Wsign-compare
endif ()

add_executable(tests bitvector_tests.cpp ../include/elias_fano.hpp eliasfano_tests.cpp rmq_tests.cpp ../include/nlogn_rmq.hpp ../include/cartesian_tree_rmq.hpp ../include/linear_space_rmq.hpp)

#target_link_libraries(tests bitvector)

target_link_libraries(tests gtest gtest_main)
