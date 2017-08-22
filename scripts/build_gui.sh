#!/usr/bin/env bash

# Goto sources directory
cd ..

# Cleanup
rm -rf build_gui

mkdir -p build_gui && cd build_gui
#../scripts/cmake.ogs.sh -DOGS_USE_QT=ON -DOGS_PACKAGING=ON ..
cmake -DOGS_USE_QT=ON -DOGS_PACKAGING=ON ..
make -j
#../scripts/cmake.ogs.sh ..
cmake ..
make -j
make package
