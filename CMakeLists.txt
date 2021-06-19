cmake_minimum_required(VERSION 3.4...3.18)

message("PYTHON VERSION=${VERSION}")

project(PyMieSim/Scatterer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python ${VERSION} COMPONENTS Interpreter Development EXACT)

add_subdirectory(extern/pybind11)

FIND_PACKAGE(PythonInterp ${VERSION} EXACT)

set(complex_bessel_DIR extern/complex_bessel)

find_package(complex_bessel)

INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIRS} )

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/PyMieSim/Tools/)
pybind11_add_module(_Coupling MODULE PyMieSim/Tools/Coupling.cpp )


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/PyMieSim/LMT)
pybind11_add_module(Scatterer MODULE PyMieSim/LMT/cpp/interface.cpp )
target_link_libraries(Scatterer PUBLIC complex_bessel::complex_bessel)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/PyMieSim/GLMT)
pybind11_add_module(_Scatterer MODULE PyMieSim/GLMT/cpp/Interface.cpp )
set_target_properties( _Scatterer PROPERTIES OUTPUT_NAME "Scatterer" )


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/PyMieSim/Tools/)
pybind11_add_module(Fibonacci MODULE PyMieSim/Tools/FibonnaciMesh.cpp )


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/PyMieSim/GLMT)
pybind11_add_module(GaussianBeam MODULE PyMieSim/GLMT/cpp/GaussianBeam.cpp )
