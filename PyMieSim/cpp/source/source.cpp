// #include "definitions.cpp" // Ensure this file provides the necessary definitions like PI and complex128
#include <vector>
#include <complex>
#include <cmath> // For std::isnan
#include "source/source.h"


double Gaussian::compute_amplitude_from_power(double wavelength, double NA, double optical_power)
{
    double omega = 0.61 * wavelength / NA;
    double area = 3.1415926535 * pow(omega / 2, 2);
    double intensity = optical_power / area;
    return sqrt(2.0 * intensity / (C_ * EPSILON0));
}