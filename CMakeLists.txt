cmake_minimum_required(VERSION 3.25)
project(ads)

set(CMAKE_CXX_STANDARD 17)

add_executable(ads src/main.cpp)

# Don't add tests and benchmarks as subdirectories because they require external libraries
#add_subdirectory(tests)
#add_subdirectory(benchmarks)
