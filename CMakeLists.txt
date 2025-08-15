cmake_minimum_required(VERSION 3.20)
set(CMAKE_OSX_DEPLOYMENT_TARGET "15.0" CACHE STRING "Minimum OS X deployment version")
project(PyMieSim LANGUAGES Fortran CXX)

# CMake settings
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)

# Set the output directory for libraries
set(LOCAL_CXX_DIR "${PROJECT_NAME}/cpp")
set(LOCAL_BIN_DIR "${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/binary")

# Compiler and linker options
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-dummy-argument -Wno-compare-reals -Wno-do-subscript -Wno-intrinsic-shadow")

# Include directories
include_directories("${LOCAL_CXX_DIR}")

# Find dependencies
if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    include_directories("/opt/homebrew/opt/libomp/include")   # for if compiler was installed through brew
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

# Print some messages
message(STATUS "OPENMP flags: ${OpenMP_CXX_FLAGS}")
message(STATUS "PyMieSim version is: ${PYMIESIM_VERSION}")
message(STATUS "PyMieSim includes directory is: ${LOCAL_CXX_DIR}")
message(STATUS "Python version to be compiled against: ${PYBIND11_PYTHON_VERSION}")
message(STATUS "Binary will be installed in location: ${LOCAL_BIN_DIR}")


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
