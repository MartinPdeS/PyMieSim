#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>


typedef std::complex<double> complex128;

#define JJ complex128(0.0,1.0)
#define PI (double)3.14159265358979323846264338
#define EPSILON0 (double)8.854187817620389e-12
#define C (double)299792458.0
