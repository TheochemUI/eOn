#!/bin/bash

echo ""
cur_dir=${PWD}

echo ""
echo " ............... Prepare the environment..."
echo ""
source modules.sh
rm -rf build
mkdir -p build

echo ""
echo " ............... Build the GTest library..."
echo ""
cd gprdimer
git submodule init
git submodule update
cd googletest/googletest
if ! patch -R -p0 -s -f --dry-run <../../cmake_gtest.patch; then
  patch < ../../cmake_gtest.patch
fi
rm -rf build
mkdir -p build
cd build
cmake -DCMAKE_CXX_FLAGS="-std=c++11" ../
make all -j 4

export LIBRARY_PATH="$cur_dir/gprdimer/googletest/googletest/build/lib:$LIBRARY_PATH"
export CPATH="$cur_dir/gprdimer/googletest/googletest/include:$CPATH"

echo ""
echo " ............... Build the EON client..."
echo ""
cd $cur_dir
rm -rf build
mkdir build && cd build

cmake .. -DCMAKE_BUILD_TYPE=Debug -DPACKAGE_TESTS=ON -DNO_WARN=TRUE \
-DGTEST_LIBRARY=$cur_dir/gprdimer/googletest/googletest/build/lib \
-DGTEST_INCLUDE_DIR=$cur_dir/gprdimer/googletest/googletest/include \
-DGTEST_MAIN_LIBRARY=$cur_dir/gprdimer/googletest/googletest/build/lib/libgtest_main.a \

make -j 24 VERBOSE=1 && make check
#export PATH=$(pwd):$PATH
