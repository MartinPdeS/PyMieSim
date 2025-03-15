cmake_minimum_required(VERSION 3.20)
set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0" CACHE STRING "Minimum OS X deployment version")
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


# Include directories
include_directories("${PYMIESIM_CXX_DIR}")
if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    include_directories("/opt/homebrew/opt/libomp/include")   # for if compiler was installed through brew
endif()


# Find dependencies
find_package(OpenMP REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)


# Platform-specific settings for static linking
if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("MinGW detected on Windows")
    set(STATIC_LINK_OPTIONS "-static")
    add_compile_options(-fopenmp)
    add_link_options(-static -fopenmp -Wl,--whole-archive -lgomp -Wl,--no-whole-archive)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")
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


function(add_pymiesim_object_target name)
    # Parse keyword arguments:
    # - SOURCES: list of source file paths (relative to ${PYMIESIM_CXX_DIR})
    # - TARGETS: list of additional targets that this object library depends on
    cmake_parse_arguments(ADD_PYMIESIM_OBJECT
        ""
        "SOURCES;TARGETS"
        ""
        ${ARGN})

    # Build a list of full paths for all source files.
    set(ALL_SOURCES "")
    foreach(src IN LISTS ADD_PYMIESIM_OBJECT_SOURCES)
        list(APPEND ALL_SOURCES "${PYMIESIM_CXX_DIR}/${src}")
    endforeach()

    # If no sources were provided, create a dummy source file.
    if (ALL_SOURCES STREQUAL "")
        set(dummy_source "${CMAKE_CURRENT_BINARY_DIR}/${name}_dummy.cpp")
        file(WRITE ${dummy_source} "// Dummy source file for object library ${name}\n")
        set(ALL_SOURCES ${dummy_source})
    endif()

    # Create the object library.
    add_library(${name} OBJECT ${ALL_SOURCES})

    # Add dependencies on the provided targets.
    foreach(dep IN LISTS ADD_PYMIESIM_OBJECT_TARGETS)
        add_dependencies(${name} ${dep})
    endforeach()

    target_link_libraries(${name} PRIVATE pybind11::module)
endfunction()


function(add_pymiesim_module name)
    # Parse keyword arguments:
    # - SOURCES: list of source file paths (relative to ${PYMIESIM_CXX_DIR})
    # - TARGETS: list of additional CMake targets to link with this module
    cmake_parse_arguments(ADD_PYMIESIM_MODULE
        ""
        "SOURCES;TARGETS"
        ""
        ${ARGN})

    # Build a list of full paths for all source files.
    set(ALL_SOURCES "")
    foreach(src IN LISTS ADD_PYMIESIM_MODULE_SOURCES)
        list(APPEND ALL_SOURCES "${PYMIESIM_CXX_DIR}/${src}")
    endforeach()

    # Create the module from the collected sources.
    pybind11_add_module(${name} MODULE ${ALL_SOURCES})
    set_target_properties(${name} PROPERTIES OUTPUT_NAME ${name})

    # Link libraries: default ones plus any extra targets provided.
    target_link_libraries(${name} PRIVATE ZBessel OpenMP::OpenMP_CXX ${ADD_PYMIESIM_MODULE_TARGETS})
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        target_link_options(${name} PRIVATE -lomp)
    endif()
endfunction()

# Build ZBessel library
add_library(ZBessel STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/libraries/amos_iso_c_fortran_wrapper.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/libraries/machine.for"
    "${CMAKE_CURRENT_SOURCE_DIR}/libraries/zbesh.for"
)
target_link_options(ZBessel PUBLIC ${STATIC_LINK_OPTIONS})

add_pymiesim_object_target(_ZBessel SOURCES "utils/bessel_subroutine.cpp")

add_pymiesim_object_target(CppProperties SOURCES "properties/properties.cpp")
add_pymiesim_object_target(CppUtils SOURCES "utils/utils.cpp")
add_pymiesim_object_target(CppVSH SOURCES "utils/VSH.cpp")
add_pymiesim_object_target(CppFibonacci SOURCES "fibonacci/fibonacci.cpp" TARGETS CppUtils)
add_pymiesim_object_target(CppFullMesh SOURCES "full_mesh/full_mesh.cpp" TARGETS CppUtils)
add_pymiesim_object_target(CppSource SOURCES "source/source.cpp" TARGETS CppUtils)
add_pymiesim_object_target(CppSphere SOURCES "sphere/sphere.cpp" TARGETS CppUtils CppVSH)
add_pymiesim_object_target(CppCoreShell SOURCES "coreshell/coreshell.cpp" TARGETS CppUtils CppVSH)
add_pymiesim_object_target(CppCylinder SOURCES "cylinder/cylinder.cpp" TARGETS CppUtils CppVSH)

add_pymiesim_object_target(CppDetector SOURCES "detector/detector.cpp" TARGETS CppUtils _ZBessel CppSphere CppCoreShell CppCylinder)

add_pymiesim_object_target(CppBaseSet SOURCES "sets/base.cpp" TARGETS CppProperties)
add_pymiesim_object_target(CppDetectorSet SOURCES "sets/detector.cpp" TARGETS CppBaseSet)
add_pymiesim_object_target(CppScattererSet SOURCES "sets/scatterer.cpp" TARGETS CppBaseSet)
add_pymiesim_object_target(CppSourceSet SOURCES "sets/source.cpp" TARGETS CppBaseSet)


# Add PyMieSim modules
add_pymiesim_module(PropertiesInterface SOURCES "properties/interface.cpp" TARGETS CppProperties)
add_pymiesim_module(FibonacciInterface SOURCES "fibonacci/interface.cpp" TARGETS CppFibonacci)
add_pymiesim_module(SourceInterface SOURCES "source/interface.cpp" TARGETS CppSource)
add_pymiesim_module(DetectorInterface SOURCES "detector/interface.cpp" TARGETS CppDetector)
add_pymiesim_module(SphereInterface SOURCES "sphere/interface.cpp" TARGETS CppSphere)
add_pymiesim_module(CylinderInterface SOURCES "cylinder/interface.cpp" TARGETS CppCylinder)
add_pymiesim_module(CoreShellInterface SOURCES "coreshell/interface.cpp" TARGETS CppCoreShell)
add_pymiesim_module(ExperimentInterface SOURCES "experiment/interface.cpp")


# add_pymiesim_module(
#     DetectorInterface
#     TARGETS CppSphere
#     SOURCES "single/includes/detectors.cpp"
#             "single/includes/sphere.cpp"
#             "single/includes/coreshell.cpp"
#             "single/includes/cylinder.cpp"
# )

# add_pymiesim_module(SphereInterface "single/interface/sphere.cpp")
# add_pymiesim_module(CylinderInterface "single/interface/cylinder.cpp")
# add_pymiesim_module(CoreShellInterface "single/interface/coreshell.cpp")
# add_pymiesim_module(ModeField "single/interface/mode_field.cpp")
# add_pymiesim_module(SourceInterface "single/interface/source.cpp")
# add_pymiesim_module(Fibonacci "single/interface/mesh.cpp")
# add_pymiesim_module(SetsInterface "experiment/interface/sets.cpp")
# add_pymiesim_module(ScattererPropertiesInterface "experiment/interface/scatterer_properties.cpp")
# add_pymiesim_module(Experiment "experiment/interface/experiment.cpp")
