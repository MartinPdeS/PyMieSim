cmake_minimum_required(VERSION 3.20)

set(CMAKE_OSX_DEPLOYMENT_TARGET "15.3" CACHE STRING "Minimum OS X deployment version")

project(PyMieSim LANGUAGES Fortran CXX)

# CMake settings
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)

# Directories
set(PYMIESIM_CXX_DIR "PyMieSim/cpp")
set(PYMIESIM_BIN_DIR "${CMAKE_SOURCE_DIR}/PyMieSim/binary")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PYMIESIM_BIN_DIR}")

# Compiler and linker options
add_compile_options(-fPIC -Wall -Wextra -Wno-uninitialized)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-dummy-argument -Wno-compare-reals -Wno-do-subscript -Wno-intrinsic-shadow")

# Find packages
find_package(OpenMP REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Include directories
include_directories("${PYMIESIM_CXX_DIR}/includes" "${PYMIESIM_CXX_DIR}/headers")

# Platform-specific settings for static linking
if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("MinGW detected on Windows")
    set(STATIC_LINK_OPTIONS "-static")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message("AppleClang compiler detected")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("GNU compiler detected")
endif()

# Print some messages
message(STATUS "OPENMP flags: ${OpenMP_CXX_FLAGS}")
message(STATUS "PyMieSim version is: ${PYMIESIM_VERSION}")
message(STATUS "PyMieSim includes directory is: ${PYMIESIM_CXX_DIR}")
message(STATUS "Python version to be compiled against: ${PYBIND11_PYTHON_VERSION}")
message(STATUS "Binary will be installed in location: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

# Build ZBessel library
add_library(ZBessel STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/libraries/amos_iso_c_fortran_wrapper.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/libraries/machine.for"
    "${CMAKE_CURRENT_SOURCE_DIR}/libraries/zbesh.for"
)
target_link_options(ZBessel PUBLIC ${STATIC_LINK_OPTIONS})

# Function to add pybind11 modules
function(add_pymiesim_module name source)
    pybind11_add_module(${name} MODULE "${PYMIESIM_CXX_DIR}/${source}")
    set_target_properties(${name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${name} PRIVATE ZBessel)
    target_link_options(${name} PUBLIC ${STATIC_LINK_OPTIONS})
endfunction()

# Add modules
add_pymiesim_module(Fibonacci "mesh_interface.cpp")
add_pymiesim_module(SourceInterface "source_interface.cpp")
add_pymiesim_module(DetectorInterface "detector_interface.cpp")
add_pymiesim_module(SphereInterface "sphere_interface.cpp")
add_pymiesim_module(CylinderInterface "cylinder_interface.cpp")
add_pymiesim_module(CoreShellInterface "coreshell_interface.cpp")
add_pymiesim_module(ModeField "mode_field_interface.cpp")
add_pymiesim_module(SetsInterface "sets_interface.cpp")
add_pymiesim_module(Experiment "experiment_interface.cpp")

# Special handling for Experiment module
if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(Experiment PRIVATE OpenMP::OpenMP_CXX)
    add_compile_options(-fopenmp)
    add_link_options(-static -fopenmp -Wl,--whole-archive -lgomp -Wl,--no-whole-archive)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(Experiment PRIVATE OpenMP::OpenMP_CXX)
else()
    target_link_libraries(Experiment PRIVATE "-Xpreprocessor -fopenmp")
endif()


# Install targets
install(
    TARGETS Fibonacci SphereInterface CylinderInterface CoreShellInterface DetectorInterface Experiment SetsInterface
    DESTINATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)
