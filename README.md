## Build

Simple execute the ./buildscript.sh script to build only the src/main.cpp file. This should generate a build/ folder
within src/, in which the resulting executable can be found. It is called 'project'.
benchmarks/run_benchmarks.sh can be used to run benchmarks.
Alternatively, simply create a build/ folder, cd into it, and call cmake with the top-level CMakeLists.txt to build
everything,
which should create the executable in build/Release/project.
This also downloads googletest and google benchmark.
