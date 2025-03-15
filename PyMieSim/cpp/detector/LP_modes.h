#pragma once

#include "utils/numpy_interface.cpp"
#include <vector>
#include <cmath>
#include "utils/bessel_subroutine.h"

#define PI (double)3.14159265358979323846264338
typedef std::complex<double> complex128;


std::vector<complex128> get_LP_mode_field(
   std::vector<double> x_coords,
   std::vector<double> y_coords,
   int azimuthal_number,
   int radial_number);

pybind11::array_t<complex128> get_LP_mode_field_py(
   std::vector<double> x_coords,
   std::vector<double> y_coords,
   int azimuthal_number,
   int radial_number);


// -