#!/bin/bash
g++ -ffast-math -O3 -mtune=native -fno-strict-aliasing -shared -fPIC -std=c++14 ./module/module.cpp ./module/utl/LambertW.cc -I./module `python3-config --includes --cflags --ldflags` -ltbb -o abcpp/dvtraitsim_cpp`python3-config --extension-suffix`
