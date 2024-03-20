cmake_minimum_required(VERSION 3.20)

file(STRINGS PyMieSim/VERSION PYMIESIM_VERSION)

project(PyMieSim LANGUAGES Fortran CXX)

message("CMAKE_CURRENT_LIST_DIR is ${CMAKE_CURRENT_LIST_DIR}")
message("CMAKE_CURRENT_SOURCE_DIR is ${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(PYMIESIM_CXX_DIR "PyMieSim/cpp")
set(PYMIESIM_BIN_DIR "${CMAKE_SOURCE_DIR}/PyMieSim/binary")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PYMIESIM_BIN_DIR}")


set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)

#FIND_PACKAGE------------------------
find_package(OpenMP)
find_package(Python ${PYBIND11_PYTHON_VERSION} EXACT COMPONENTS Interpreter Development)
add_subdirectory("extern/pybind11")

include_directories("${CMAKE_SOURCE_DIR}/PyMieSim/cpp/includes")

add_compile_options(-Wno-return-type)

if (WIN32)
    message("Windows OS detected")
    set(STATIC_STDLIB "-w -lm -static-libgcc -static-libstdc++ -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive") # Need to keep -lwinpthread static because of MinGW dll directory not the same as windows' MSVC
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message("AppleClang compiler detected")
    set(STATIC_STDLIB "-w -lm -static-libstdc++")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("GNU compiler detected")
    set(STATIC_STDLIB "-w -lm -static-libstdc++ -static-libgcc")
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
target_link_libraries(Fibonacci PRIVATE "${STATIC_STDLIB}")

#Building SPHERE-----------------------------------------------------
pybind11_add_module(CXX_SPHERE MODULE "${PYMIESIM_CXX_DIR}/sphere_interface.cpp" )
set_target_properties(CXX_SPHERE PROPERTIES OUTPUT_NAME "SphereInterface")
target_link_libraries(CXX_SPHERE PRIVATE ZBessel "${STATIC_STDLIB}")

#Building CYLINDER---------------------------------------------------
pybind11_add_module(CXX_CYLINDER MODULE "${PYMIESIM_CXX_DIR}/cylinder_interface.cpp" )
set_target_properties(CXX_CYLINDER PROPERTIES OUTPUT_NAME "CylinderInterface")
target_link_libraries(CXX_CYLINDER PRIVATE ZBessel "${STATIC_STDLIB}")

#Building CORESHELL--------------------------------------------------
pybind11_add_module(CXX_CORESHELL MODULE "${PYMIESIM_CXX_DIR}/coreshell_interface.cpp" )
set_target_properties(CXX_CORESHELL PROPERTIES OUTPUT_NAME "CoreShellInterface")
target_link_libraries(CXX_CORESHELL PRIVATE ZBessel "${STATIC_STDLIB}")

#Building DETECTOR--------------------------------------------------
pybind11_add_module(CXX_DETECTOR MODULE "${PYMIESIM_CXX_DIR}/detector_interface.cpp" )
set_target_properties(CXX_DETECTOR PROPERTIES OUTPUT_NAME "DetectorInterface")
target_link_libraries(CXX_DETECTOR PRIVATE "${STATIC_STDLIB}")

#Building EXPERIMENT------------------------------------------------
if (WIN32)
    pybind11_add_module(Experiment MODULE "${PYMIESIM_CXX_DIR}/experiment_interface.cpp" )
    set_target_properties(Experiment PROPERTIES OUTPUT_NAME "Experiment")
    target_link_libraries(Experiment PRIVATE ZBessel "${STATIC_STDLIB}")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    pybind11_add_module(Experiment MODULE "${PYMIESIM_CXX_DIR}/experiment_interface.cpp" )
    set_target_properties(Experiment PROPERTIES OUTPUT_NAME "Experiment")
    target_link_libraries(Experiment PRIVATE ZBessel "${STATIC_STDLIB}" OpenMP::OpenMP_CXX )
else(APPLE)
    pybind11_add_module(Experiment MODULE "${PYMIESIM_CXX_DIR}/experiment_interface.cpp" )
    set_target_properties(Experiment PROPERTIES OUTPUT_NAME "Experiment")
    target_link_libraries(Experiment PRIVATE ZBessel "${STATIC_STDLIB} -Xpreprocessor -fopenmp")
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


