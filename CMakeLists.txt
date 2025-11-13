cmake_minimum_required(VERSION 3.20)
set(CMAKE_OSX_DEPLOYMENT_TARGET "15.0" CACHE STRING "Minimum OS X deployment version")
project(PyMieSim LANGUAGES Fortran CXX)

# CMake settings
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build")

# Set the output directory for libraries
set(LOCAL_CXX_DIR "${PROJECT_NAME}/cpp")
set(LOCAL_BIN_DIR "${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/binary")

# Compiler and linker options
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-dummy-argument -Wno-compare-reals -Wno-do-subscript -Wno-intrinsic-shadow -Wmaybe-uninitialized")

# Include directories
include_directories("${LOCAL_CXX_DIR}")

# Find dependencies
if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "")
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "")
    find_library(OpenMP_omp_LIBRARY omp HINTS $ENV{LIBOMP_PREFIX}/lib /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
    if(OpenMP_omp_LIBRARY)
        set(OpenMP_CXX_LIBRARIES "${OpenMP_omp_LIBRARY}" CACHE STRING "")
    endif()
endif()

find_package(OpenMP REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Platform-specific settings for static linking
if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("MinGW detected on Windows")
    set(STATIC_LINK_OPTIONS "-static")
    add_compile_options(-fopenmp)
    add_link_options(-static -fopenmp -Wl,--whole-archive -lgomp -Wl,--no-whole-archive)
endif()

# Compiler and linker options
add_compile_options(-Wall -Wextra -pedantic-errors)



message(STATUS "")
message(STATUS "==================== Build configuration ====================")

message(STATUS "Compiler information")
message(STATUS "  C compiler           : ${CMAKE_C_COMPILER}")
message(STATUS "  C compiler ID        : ${CMAKE_C_COMPILER_ID}")
message(STATUS "  C compiler version   : ${CMAKE_C_COMPILER_VERSION}")

message(STATUS "  CXX compiler         : ${CMAKE_CXX_COMPILER}")
message(STATUS "  CXX compiler ID      : ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "  CXX compiler version : ${CMAKE_CXX_COMPILER_VERSION}")

message(STATUS "  Fortran compiler     : ${CMAKE_Fortran_COMPILER}")
message(STATUS "  Fortran ID           : ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "  Fortran version      : ${CMAKE_Fortran_COMPILER_VERSION}")

message(STATUS "")
message(STATUS "OpenMP configuration")
message(STATUS "  OpenMP flags         : ${OpenMP_CXX_FLAGS}")
message(STATUS "  OpenMP detected      : ${OpenMP_FOUND}")
message(STATUS "  OpenMP version       : ${OpenMP_CXX_VERSION}")

message(STATUS "")
message(STATUS "Python configuration")
message(STATUS "  Python executable    : ${Python_EXECUTABLE}")
message(STATUS "  Python version       : ${PYBIND11_PYTHON_VERSION}")
message(STATUS "  Python include dir   : ${Python_INCLUDE_DIRS}")
message(STATUS "  Python site packages : ${Python_SITELIB}")

message(STATUS "")
message(STATUS "PyMieSim integration")
message(STATUS "  PyMieSim version     : ${PYMIESIM_VERSION}")
message(STATUS "  PyMieSim include dir : ${LOCAL_CXX_DIR}")

message(STATUS "")
message(STATUS "Install destination")
message(STATUS "  Binary output path   : ${LOCAL_BIN_DIR}")

message(STATUS "==============================================================")
message(STATUS "")


add_subdirectory(PyMieSim/cpp/coordinates)                  # coordinates
add_subdirectory(PyMieSim/cpp/bessel_subroutine)            # bessel_subroutine
add_subdirectory(PyMieSim/cpp/mode_field)                   # cpp_mode_field
add_subdirectory(PyMieSim/cpp/full_mesh)                    # cpp_full_mesh
add_subdirectory(PyMieSim/cpp/fibonacci)                    # cpp_fibonacci      - interface_fibonacci
add_subdirectory(PyMieSim/cpp/source)                       # cpp_source         - interface_source
add_subdirectory(PyMieSim/cpp/scatterer)                    # cpp_scatterers
add_subdirectory(PyMieSim/cpp/detector)                     # cpp_detector       - interface_detector
add_subdirectory(PyMieSim/cpp/sets)                         # cpp_sets           - interface_sets
add_subdirectory(PyMieSim/cpp/experiment)                   # cpp_experiment     - interface_experiment
