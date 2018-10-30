#!/bin/bash
module load tbb
module load Python/3.6.4-foss-2018a
g++ -ffast-math -O3 -march=znver1 -fno-strict-aliasing -shared -fPIC -std=c++14 ./module/module.cpp ./module/utl/LambertW.cc -I./module `python3-config --includes --cflags --ldflags` -ltbb -o abcpp/dvtraitsim_cpp`python3-config --extension-suffix`
