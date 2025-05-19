#!/bin/bash

# Source and output
SRC="Calc_CPP_Profiler.cpp Conc_Calc_CPP.cpp Params.cpp Reac_Rate_CPP.cpp"
OUT="calc_cpp_app_profiler"

# Include paths
INCLUDE_DIRS="-I/opt/homebrew/opt/gperftools/include"

# Library paths and libraries
LIB_DIRS="-L/opt/homebrew/opt/gperftools/lib"
LIBS="-lprofiler"

# Compile using clang++
clang++ -std=c++11 -Wall -g -D PROFILING $SRC "$INCLUDE_DIRS" "$LIB_DIRS" "$LIBS" -o $OUT #-g keeps variable names

# Check compilation result
if [ $? -eq 0 ]; then
    echo "Compilation successful. Output: $OUT"
else
    echo "Compilation failed."
    exit 1
fi

codesign -s - -v -f --entitlements /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/debug.plist /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/calc_cpp_app_profiler

codesign -dvvv --entitlements - /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/calc_cpp_app_profiler

