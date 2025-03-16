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


function(add_pymiesim_library_target name)
    # Parse keyword arguments:
    # - SOURCES: list of source file paths (relative to ${PYMIESIM_CXX_DIR})
    # - TARGETS: list of additional targets that this library depends on
    # - TYPE: optional type specifier, either SHARED or STATIC (default: STATIC)
    cmake_parse_arguments(ADD_PYMIESIM_LIB
        "SHARED"
        "SOURCES;TARGETS;TYPE"
        ""
        ${ARGN})

    # Build a list of full paths for all source files.
    set(ALL_SOURCES "")
    foreach(src IN LISTS ADD_PYMIESIM_LIB_SOURCES)
        list(APPEND ALL_SOURCES "${PYMIESIM_CXX_DIR}/${src}")
    endforeach()

    # If no sources were provided, create a dummy source file.
    if (ALL_SOURCES STREQUAL "")
        set(dummy_source "${CMAKE_CURRENT_BINARY_DIR}/${name}_dummy.cpp")
        file(WRITE ${dummy_source} "// Dummy source file for library ${name}\n")
        set(ALL_SOURCES ${dummy_source})
    endif()

    # Determine library type: If the TYPE argument is provided, use that.
    # Otherwise, use SHARED if the SHARED keyword was set, else default to STATIC.
    if (ADD_PYMIESIM_LIB_TYPE)
        set(LIB_TYPE ${ADD_PYMIESIM_LIB_TYPE})
    elseif(ADD_PYMIESIM_LIB_SHARED)
        set(LIB_TYPE SHARED)
    else()
        set(LIB_TYPE STATIC)
    endif()

    # Create the library.
    add_library(${name} ${LIB_TYPE} ${ALL_SOURCES})

    # Add dependencies on the provided targets.
    foreach(dep IN LISTS ADD_PYMIESIM_LIB_TARGETS)
        target_link_libraries(${name} PRIVATE ${dep} pybind11::module)
    endforeach()

    # Optionally, add any common include directories or compile definitions here.
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
add_pymiesim_library_target(_ZBessel SOURCES "utils/bessel_subroutine.cpp" TYPE SHARED)

add_subdirectory(PyMieSim/cpp/mode_field)
add_subdirectory(PyMieSim/cpp/full_mesh)
add_subdirectory(PyMieSim/cpp/fibonacci)
add_subdirectory(PyMieSim/cpp/source)
add_subdirectory(PyMieSim/cpp/sphere)
add_subdirectory(PyMieSim/cpp/cylinder)
add_subdirectory(PyMieSim/cpp/coreshell)
add_subdirectory(PyMieSim/cpp/detector)

# add_pymiesim_library_target(_ZBessel SOURCES "utils/bessel_subroutine.cpp" TYPE SHARED)
# add_pymiesim_library_target(CppProperties SOURCES "properties/properties.cpp" TYPE SHARED)
# add_pymiesim_library_target(CppUtils SOURCES "utils/utils.cpp" TYPE SHARED)
# add_pymiesim_library_target(CppVSH SOURCES "utils/VSH.cpp" TYPE SHARED)
# add_pymiesim_library_target(CppFibonacci SOURCES "fibonacci/fibonacci.cpp" TARGETS CppUtils TYPE SHARED)
# add_pymiesim_library_target(CppFullMesh SOURCES "full_mesh/full_mesh.cpp" TARGETS CppUtils TYPE SHARED)
# add_pymiesim_library_target(CppSource SOURCES "source/source.cpp" TARGETS CppUtils TYPE SHARED)
# add_pymiesim_library_target(CppSphere SOURCES "sphere/sphere.cpp" TARGETS CppUtils CppVSH TYPE SHARED)
# add_pymiesim_library_target(CppCoreShell SOURCES "coreshell/coreshell.cpp" TARGETS CppUtils CppVSH TYPE SHARED)
# add_pymiesim_library_target(CppCylinder SOURCES "cylinder/cylinder.cpp" TARGETS CppUtils CppVSH TYPE SHARED)
# add_pymiesim_library_target(CppDetector SOURCES "detector/detector.cpp" TARGETS CppUtils _ZBessel CppSphere CppCoreShell CppCylinder TYPE SHARED)
# add_pymiesim_library_target(CppBaseSet SOURCES "sets/base.cpp" TARGETS CppProperties TYPE SHARED)
# add_pymiesim_library_target(CppDetectorSet SOURCES "sets/detector.cpp" TARGETS CppBaseSet TYPE SHARED)
# add_pymiesim_library_target(CppScattererSet SOURCES "sets/scatterer.cpp" TARGETS CppBaseSet TYPE SHARED)
# add_pymiesim_library_target(CppSourceSet SOURCES "sets/source.cpp" TARGETS CppBaseSet TYPE SHARED)


# Add PyMieSim modules
# add_pymiesim_module(PropertiesInterface SOURCES "properties/interface.cpp" TARGETS CppProperties)
# add_pymiesim_module(FibonacciInterface SOURCES "fibonacci/interface.cpp" TARGETS CppFibonacci)
# add_pymiesim_module(SourceInterface SOURCES "source/interface.cpp" TARGETS CppSource)
# add_pymiesim_module(DetectorInterface SOURCES "detector/interface.cpp" TARGETS CppDetector CppUtils _ZBessel CppSphere CppCoreShell CppCylinder)
# add_pymiesim_module(SphereInterface SOURCES "sphere/interface.cpp" TARGETS CppSphere)
# add_pymiesim_module(CylinderInterface SOURCES "cylinder/interface.cpp" TARGETS CppCylinder)
# add_pymiesim_module(CoreShellInterface SOURCES "coreshell/interface.cpp" TARGETS CppCoreShell)
# add_pymiesim_module(ExperimentInterface SOURCES "experiment/interface.cpp")
