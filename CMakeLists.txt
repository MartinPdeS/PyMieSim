cmake_minimum_required(VERSION 3.20)
set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0" CACHE STRING "Minimum OS X deployment version")
project(PyMieSim LANGUAGES Fortran CXX)

# CMake settings
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)


# Set the output directory for libraries
set(PYMIESIM_CXX_DIR "PyMieSim/cpp")
set(PYMIESIM_BIN_DIR "${CMAKE_SOURCE_DIR}/PyMieSim/binary")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PYMIESIM_BIN_DIR}")


# Compiler and linker options
add_compile_options(-fPIC -Wall -Wextra -Wno-uninitialized)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-dummy-argument -Wno-compare-reals -Wno-do-subscript -Wno-intrinsic-shadow")


# Include directories
include_directories("${PYMIESIM_CXX_DIR}")
if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    include_directories("/opt/homebrew/opt/libomp/include")   # for if compiler was installed through brew
endif()


# Find dependencies
find_package(OpenMP REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)


# Print some messages
message(STATUS "OPENMP flags: ${OpenMP_CXX_FLAGS}")
message(STATUS "PyMieSim version is: ${PYMIESIM_VERSION}")
message(STATUS "PyMieSim includes directory is: ${PYMIESIM_CXX_DIR}")
message(STATUS "Python version to be compiled against: ${PYBIND11_PYTHON_VERSION}")
message(STATUS "Binary will be installed in location: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")


# add_subdirectory(PyMieSim/cpp/_bessel)
add_subdirectory(PyMieSim/cpp/bessel)
add_subdirectory(PyMieSim/cpp/mode_field)         # cpp_mode_field
add_subdirectory(PyMieSim/cpp/full_mesh)          # cpp_full_mesh
add_subdirectory(PyMieSim/cpp/fibonacci)
add_subdirectory(PyMieSim/cpp/source)             # cpp_source
add_subdirectory(PyMieSim/cpp/sphere)
add_subdirectory(PyMieSim/cpp/cylinder)
add_subdirectory(PyMieSim/cpp/coreshell)
add_subdirectory(PyMieSim/cpp/detector)           # cpp_detector
add_subdirectory(PyMieSim/cpp/sets)
add_subdirectory(PyMieSim/cpp/experiment)         # cpp_experiment


set(ALL_TARGETS cpp_source cpp_full_mesh cpp_fibonacci cpp_mode_field
                cpp_sphere cpp_cylinder cpp_coreshell cpp_detector cpp_sets cpp_experiment
                interface_source interface_mode_field interface_fibonacci interface_sphere
                interface_cylinder interface_coreshell interface_detector interface_sets interface_experiment)

install(TARGETS ${ALL_TARGETS}
        LIBRARY DESTINATION ${PYMIESIM_BIN_DIR}
        ARCHIVE DESTINATION ${PYMIESIM_BIN_DIR}
        RUNTIME DESTINATION ${PYMIESIM_BIN_DIR})