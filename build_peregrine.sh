#!/bin/bash
module load Python/3.6.4-intel-2018a  # intel's python3 
module load tbb/2018_U5-GCCcore-6.4.0
icc -O3 -shared -fno-strict-aliasing -fPIC -std=c++14 -DNDEBUG ./module/module.cpp ./module/utl/LambertW.cc -I./module `python3-config --includes --cflags --ldflags` -ltbb -ltbbmalloc -o abcpp/dvtraitsim_cpp`python3-config --extension-suffix`
