cmake_minimum_required(VERSION 3.25)

project(project)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

if (MSVC)
    #    TODO
else ()
    add_compile_options(-Ofast -march=native -flto=auto -fno-rtti -fno-exceptions -fvisibility=hidden)
    add_link_options(-flto=auto)
endif ()

# uncomment if the machine has inefficient support for BMI2 intructions like pdep, like AMD processors before Zen 3.
add_compile_options(-DADS_USE_BMI2_INTRINSICS)
add_compile_options(-DADS_NO_FALLBACKS)

add_executable(project main.cpp)
