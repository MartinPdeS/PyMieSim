#pragma once

#include <vector>
#include <complex>
#include <cmath> // For std::isnan

#define PI (double)3.14159265358979323846264338
#define EPSILON0 (double)8.854187817620389e-12
#define C_ (double)299792458.0
typedef std::complex<double> complex128;


class BaseSource {
    public:
        double wavelength;
        std::vector<complex128> jones_vector;
        double amplitude;
        double wavenumber;
        std::vector<size_t> indices;
        size_t wavelength_index;

    BaseSource() = default;
    BaseSource(double wavelength, std::vector<complex128> jones_vector, double amplitude);

};

class Planewave: public BaseSource {
    public:
        Planewave() = default;

        Planewave(double wavelength, std::vector<complex128> jones_vector, double amplitude);
    };

class Gaussian: public BaseSource {
    public:
        double NA;
        double optical_power;

        Gaussian() = default;

        Gaussian(double wavelength, std::vector<complex128> jones_vector, double NA, double optical_power);

        double compute_amplitude_from_power(double wavelength, double NA, double optical_power);
    };

