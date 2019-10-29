#!/bin/bash
module load Python/3.6.4-foss-2019a
module load tbb/2018_U5-GCCcore-6.4.0
g++ -O3 -march=znver1 -mavx2 -fno-strict-aliasing -shared -fPIC -fexceptions -std=c++14 ./module/module.cpp ./module/utl/LambertW.cc -I./module `python3-config --includes --cflags --ldflags` -ltbb -ltbbmalloc -o abcpp/dvtraitsim_cpp`python3-config --extension-suffix`
