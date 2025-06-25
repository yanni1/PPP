#!/bin/bash

#source and output
SRC="Calc_CPP_Profiler.cpp Conc_Calc_CPP.cpp Params.cpp Reac_Rate_CPP.cpp"
OUT="calc_cpp_app_profiler"

#compile with clang++
clang++ -std=c++11 -Wall -g -D PROFILING $SRC -o $OUT #-g keeps variable names

#check compilation result
if [ $? -eq 0 ]; then
    echo "Compilation successful. Output: $OUT"
else
    echo "Compilation failed."
    exit 1
fi

codesign -s - -v -f --entitlements /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/debug.plist /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/calc_cpp_app_profiler

codesign -dvvv --entitlements - /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/calc_cpp_app_profiler

