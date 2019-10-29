#!/bin/bash
g++ -O3 -mtune=native -mavx2 -fno-strict-aliasing -shared -fexceptions -fPIC -std=c++14 -DNDEBUG ./module/module.cpp ./module/utl/LambertW.cc -I./module `python3-config --includes --cflags --ldflags` -ltbb -ltbbmalloc -o abcpp/dvtraitsim_cpp`python3-config --extension-suffix`
