#pragma once

#include <vector>
#include <complex>
#include "utils/numpy_interface.h"
#include <pybind11/numpy.h>

#define PI (double)3.14159265358979323846264338
typedef std::complex<double> complex128;



// ________________________________ LP ___________________________________

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

// ________________________________ HG ___________________________________


double hermite_next(unsigned n, double x, double Hn, double Hnm1);

double hermite_imp(unsigned n, double x);

// Helper function to calculate the Hermite-Gaussian mode field amplitude
std::vector<complex128> get_HG_mode_field(
   std::vector<double> x_coords,
   std::vector<double> y_coords,
   int x_number,
   int y_number,
   double wavelength = 1.55,
   double waist_radius = 0.3,
   double z = 0.0);

// Helper function to calculate the Hermite-Gaussian mode field amplitude
pybind11::array_t<complex128> get_HG_mode_field_py(
   std::vector<double> x_coords,
   std::vector<double> y_coords,
   int x_number,
   int y_number);

// _____________________________ LG__________________________________


double laguerre_next(unsigned n, double x, double Ln, double Lnm1);


double laguerre_next(unsigned n, unsigned l, double x, double Pl, double Plm1);


double laguerre_imp(unsigned n, double x);


double laguerre_imp(unsigned n, unsigned m, double x);


std::vector<complex128> get_LG_mode_field(
   std::vector<double> x_coords,
   std::vector<double> y_coords,
   int azimuthal_number,
   int radial_number,
   double wavelength = 1.55,
   double waist_radius = 0.3,
   double z = 0.0);

pybind11::array_t<complex128> get_LG_mode_field_py(
   std::vector<double> x_coords,
   std::vector<double> y_coords,
   int azimuthal_number,
   int radial_number);

