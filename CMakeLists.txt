cmake_minimum_required(VERSION 3.20..)

set(CMAKE_OSX_DEPLOYMENT_TARGET "15.3" CACHE STRING "Minimum OS X deployment version")

project(PyMieSim LANGUAGES Fortran CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(PYMIESIM_CXX_DIR "PyMieSim/cpp")
set(PYMIESIM_BIN_DIR "${CMAKE_SOURCE_DIR}/PyMieSim/binary")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PYMIESIM_BIN_DIR}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-dummy-argument -Wno-compare-reals -Wno-do-subscript -Wno-intrinsic-shadow")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)

# Ensure -fPIC is used globally
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_options(-fPIC)

#FIND_PACKAGE------------------------
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(OpenMP)

include_directories("${CMAKE_SOURCE_DIR}/PyMieSim/cpp/includes")
include_directories("${CMAKE_SOURCE_DIR}/PyMieSim/cpp/headers")

add_compile_options(-Wall -Wextra -Wno-uninitialized)

# Platform-specific settings for static linking
if (WIN32)
    message("Windows OS detected")
    set(STATIC_LINK_OPTIONS "-static" "-static-libgcc" "-static-libstdc++")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message("AppleClang compiler detected")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("GNU compiler detected")
endif()


message(STATUS "OPENMP flags: "  ${OpenMP_CXX_FLAGS})
message(STATUS "PyMieSim  version is: " ${PYMIESIM_VERSION})
message(STATUS "PyMieSim includes directory is: " ${PYMIESIM_CXX_DIR})
message(STATUS "Python version to being compiled against: " ${PYBIND11_PYTHON_VERSION})
message(STATUS "Binary will be installed in location: " ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "Custom flags passed to linker: " ${STATIC_STDLIB})

#Building libZBessel-------------------------------------------------
add_library(
    ZBessel STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/libraries/amos_iso_c_fortran_wrapper.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/libraries/machine.for
    ${CMAKE_CURRENT_SOURCE_DIR}/libraries/zbesh.for
)
target_link_libraries(ZBessel PUBLIC "${STATIC_STDLIB}")


#Building FibonacciMesh----------------------------------------------
pybind11_add_module(Fibonacci MODULE "${PYMIESIM_CXX_DIR}/mesh_interface.cpp")
set_target_properties(Fibonacci PROPERTIES OUTPUT_NAME "Fibonacci")
target_link_libraries(Fibonacci PRIVATE ZBessel " ")
target_link_options(Fibonacci PUBLIC ${STATIC_LINK_OPTIONS})

#Building SOURCE--------------------------------------------------
pybind11_add_module(CXX_SOURCE MODULE "${PYMIESIM_CXX_DIR}/source_interface.cpp" )
set_target_properties(CXX_SOURCE PROPERTIES OUTPUT_NAME "SourceInterface")
target_link_libraries(CXX_SOURCE PRIVATE ZBessel " ")
target_link_options(CXX_SOURCE PUBLIC ${STATIC_LINK_OPTIONS})

#Building DETECTOR--------------------------------------------------
pybind11_add_module(CXX_DETECTOR MODULE "${PYMIESIM_CXX_DIR}/detector_interface.cpp" )
set_target_properties(CXX_DETECTOR PROPERTIES OUTPUT_NAME "DetectorInterface")
target_link_libraries(CXX_DETECTOR PRIVATE ZBessel " ")
target_link_options(CXX_DETECTOR PUBLIC ${STATIC_LINK_OPTIONS})

#Building SPHERE-----------------------------------------------------
pybind11_add_module(CXX_SPHERE MODULE "${PYMIESIM_CXX_DIR}/sphere_interface.cpp" )
set_target_properties(CXX_SPHERE PROPERTIES OUTPUT_NAME "SphereInterface")
target_link_libraries(CXX_SPHERE PRIVATE ZBessel " ")
target_link_options(CXX_SPHERE PUBLIC ${STATIC_LINK_OPTIONS})

#Building CYLINDER---------------------------------------------------
pybind11_add_module(CXX_CYLINDER MODULE "${PYMIESIM_CXX_DIR}/cylinder_interface.cpp" )
set_target_properties(CXX_CYLINDER PROPERTIES OUTPUT_NAME "CylinderInterface")
target_link_libraries(CXX_CYLINDER PRIVATE ZBessel " ")
target_link_options(CXX_CYLINDER PUBLIC ${STATIC_LINK_OPTIONS})

#Building CORESHELL--------------------------------------------------
pybind11_add_module(CXX_CORESHELL MODULE "${PYMIESIM_CXX_DIR}/coreshell_interface.cpp" )
set_target_properties(CXX_CORESHELL PROPERTIES OUTPUT_NAME "CoreShellInterface")
target_link_libraries(CXX_CORESHELL PRIVATE ZBessel " ")
target_link_options(CXX_CORESHELL PUBLIC ${STATIC_LINK_OPTIONS})

#Building MODEFIELD--------------------------------------------------
pybind11_add_module(CXX_MODEFIELD MODULE "${PYMIESIM_CXX_DIR}/mode_field_interface.cpp" )
set_target_properties(CXX_MODEFIELD PROPERTIES OUTPUT_NAME "ModeField")
target_link_options(CXX_MODEFIELD PUBLIC ${STATIC_LINK_OPTIONS})

#Building EXPERIMENT------------------------------------------------
if (WIN32)
    pybind11_add_module(Experiment MODULE "${PYMIESIM_CXX_DIR}/experiment_interface.cpp" )
    set_target_properties(Experiment PROPERTIES OUTPUT_NAME "Experiment")
    target_link_libraries(Experiment PRIVATE ZBessel " ")
    target_link_options(Experiment PUBLIC ${STATIC_LINK_OPTIONS})
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    pybind11_add_module(Experiment MODULE "${PYMIESIM_CXX_DIR}/experiment_interface.cpp" )
    set_target_properties(Experiment PROPERTIES OUTPUT_NAME "Experiment")
    target_link_libraries(Experiment PRIVATE ZBessel OpenMP::OpenMP_CXX )
    target_link_options(Experiment PUBLIC ${STATIC_LINK_OPTIONS})
else(APPLE)
    pybind11_add_module(Experiment MODULE "${PYMIESIM_CXX_DIR}/experiment_interface.cpp" )
    set_target_properties(Experiment PROPERTIES OUTPUT_NAME "Experiment")
    target_link_libraries(Experiment PRIVATE ZBessel "-Xpreprocessor -fopenmp")
    target_link_options(Experiment PUBLIC ${STATIC_LINK_OPTIONS})
endif()

#Building CppSets----------------------------------------------------
pybind11_add_module(Sets MODULE "${PYMIESIM_CXX_DIR}/sets_interface.cpp" )
set_target_properties(Sets PROPERTIES OUTPUT_NAME "Sets")
target_link_libraries(Sets PRIVATE ZBessel "${STATIC_STDLIB}")

set_target_properties(ZBessel Fibonacci CXX_SPHERE CXX_CYLINDER CXX_CORESHELL CXX_DETECTOR Experiment Sets
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    POSITION_INDEPENDENT_CODE TRUE
)


install(TARGETS Fibonacci CXX_SPHERE CXX_CYLINDER CXX_CORESHELL CXX_DETECTOR Experiment Sets
        # RUNTIME_DEPENDENCIES # do not uncomment this as it will mixed dependencies that may not be manylinux-compatible
        DESTINATION "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
        )


