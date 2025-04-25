project(Conc_Calc_CPP)
cmake_minimum_required(VERSION 3.18...3.22)
#set(Python_EXECUTABLE "/Users/yanni/.pyenv/versions/3.12.2/bin/python")
find_package(Python 3.12.2 COMPONENTS Interpreter Development.Module REQUIRED) # minimal python is 3.8
message(STATUS "Python found: ${Python_EXECUTABLE}")
message(STATUS "Python version: ${Python_VERSION}")
message(STATUS "Python include dir: ${Python_INCLUDE_DIRS}")
message(STATUS "Python library: ${Python_LIBRARIES}")

# Configure CMake to perform an optimized release build by default unless another build type is specified.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Detect the installed nanobind package and import it into CMake
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)
message(STATUS "Nanobind directory found: ${NB_DIR}")
message(STATUS "CMAKE_PREFIX_PATH includes: ${CMAKE_PREFIX_PATH}")

# Create the module
nanobind_add_module(Conc_Calc_CPP Conc_Calc_CPP.cpp)

# install the module
install(TARGETS Conc_Calc_CPP DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/..")

# generate/build/install _cmake_build directory
# > cmake -S . -B _cmake_build
# > cmake --build _cmake_build
# > cmake --install _cmake_build
# or
# > cmake -S . -B _cmake_build && cmake --build _cmake_build && cmake --install _cmake_build

