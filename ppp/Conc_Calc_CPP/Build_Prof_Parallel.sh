#!/bin/bash

# Source and output
SRC="Calc_CPP_Profiler.cpp Conc_Calc_CPP.cpp Params.cpp Reac_Rate_CPP.cpp"
OUT="calc_cpp_app_profiler"

# LLVM paths
LLVM_PATH="/opt/homebrew/opt/llvm"
OMP_FLAGS="-Xpreprocessor -fopenmp -I$LLVM_PATH/include -L$LLVM_PATH/lib -lomp"

# Profiler include and libs
GPERF_INCLUDE="-I/opt/homebrew/opt/gperftools/include"
GPERF_LIB="-L/opt/homebrew/opt/gperftools/lib -lprofiler"

# Add optimization and SIMD flags
OPT_FLAGS="-O3 -march=native -funroll-loops -Rpass=loop-vectorize"

# Use LLVM clang++
$LLVM_PATH/bin/clang++ -std=c++17 -Wall -g -D PROFILING $OPT_FLAGS $SRC $OMP_FLAGS $GPERF_INCLUDE $GPERF_LIB -o $OUT

# Check compilation result
if [ $? -eq 0 ]; then
    echo "Compilation successful. Output: $OUT"
else
    echo "Compilation failed."
    exit 1
fi

# Code signing
codesign -s - -v -f --entitlements /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/debug.plist /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/$OUT

codesign -dvvv --entitlements - /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/$OUT
