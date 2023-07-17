#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

cd "${SCRIPTPATH}"/src &&
mkdir -p build &&
cd build &&
cmake -DCMAKE_BUILD_TYPE=Release ../ &&
cmake --build .
