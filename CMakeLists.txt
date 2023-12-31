cmake_minimum_required(VERSION 3.25)

cmake_policy(SET CMP0069 NEW)

project(ads)

set(CMAKE_CXX_STANDARD_REQUIRED 17)

set(CMAKE_CXX_STANDARD 20)

if (MSVC)
    add_compile_options(/W4 /WX)
else ()
    add_compile_options(-Wall -Wextra -Wpedantic -march=native -fno-rtti -fno-exceptions -fvisibility=hidden)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-stdlib=libc++) # at least for C++20 features, clang doesn't fully support libstdc++
        add_link_options(-stdlib=libc++)
    endif ()
endif ()


# Don't add tests and benchmarks as subdirectories because they require external libraries?
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(benchmarks)
