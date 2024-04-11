#pragma once

#include <pybind11/numpy.h>

typedef std::complex<double> complex128;

#define JJ complex128(0.0,1.0)

#define PI (double)3.14159265358979323846264338

double EPSILON0 = 8.854187817620389e-12;
double C = 299792458.0 ;
