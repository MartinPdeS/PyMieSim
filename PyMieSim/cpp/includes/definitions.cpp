#ifndef DEFINITION_H
#define DEFINITION_H


#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

typedef unsigned int                     uint;
typedef std::complex<double>             complex128;

typedef std::vector<size_t>              IVector;
typedef std::vector<double>              DVector;
typedef std::vector<complex128>          CVector;

typedef pybind11::array_t<double>        ndarray;
typedef pybind11::array_t<complex128>    Cndarray;
typedef std::vector<std::vector<double>> Matrix3;

typedef pybind11::array_t<size_t>        Indarray;
typedef pybind11::buffer_info            info;

#define JJ complex128(0.0,1.0)

#define PI (double)3.14159265358979323846264338

double EPSILON0 = 8.854187817620389e-12;
double C = 299792458.0 ;


#endif