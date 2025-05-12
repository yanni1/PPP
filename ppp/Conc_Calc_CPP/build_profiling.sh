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

./calc_cpp_app_profiler

pprof --pdf calc_cpp_app_profiler Calc_CPP_profile.log >> /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Conc_Calc_CPP/pprof_pdf.pdf
