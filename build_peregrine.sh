#!/bin/bash
module load icc     # intel c++ compiler 
module load Python  # intel's python3 
icc -fast -O3 -shared -fno-strict-aliasing -fPIC -std=c++14 ./module/module.cpp ./module/utl/LambertW.cc -I./module `python3-config --includes --cflags --ldflags` -ltbb -o abcpp/dvtraitsim_cpp`python3-config --extension-suffix`
