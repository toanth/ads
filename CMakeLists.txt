cmake_minimum_required(VERSION 3.25)

project(ads)

set(CMAKE_CXX_STANDARD_REQUIRED 17)

set(CMAKE_CXX_STANDARD 20)

# Don't add tests and benchmarks as subdirectories because they require external libraries?
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(benchmarks)
