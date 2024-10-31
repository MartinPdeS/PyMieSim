#pragma once

#include <vector>
#include <complex>

#define PI (double)3.14159265358979323846264338
#define EPSILON0 (double)8.854187817620389e-12
#define C_ (double)299792458.0
typedef std::complex<double> complex128;

namespace SOURCE
{
    class BaseSource {
        public:
            double wavelength;
            std::vector<complex128> jones_vector;
            double amplitude;
            double wavenumber;
            std::vector<size_t> indices;
            size_t wavelength_index;

        BaseSource() = default;
        BaseSource(double wavelength, std::vector<complex128> jones_vector, double amplitude)
        : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude), wavenumber(2 * PI / wavelength)
        {}

    };

    class Planewave: public BaseSource {
        public:
            Planewave() = default;

            Planewave(double wavelength, std::vector<complex128> jones_vector, double amplitude)
            : BaseSource(wavelength, jones_vector, amplitude){}
        };

    class Gaussian: public BaseSource {
        public:
            double NA;
            double optical_power;

            Gaussian() = default;

            Gaussian(double wavelength, std::vector<complex128> jones_vector, double NA, double optical_power)
            : BaseSource(wavelength, jones_vector, compute_amplitude_from_power(wavelength, NA, optical_power)), NA(NA), optical_power(optical_power)
            {}

            static double compute_amplitude_from_power(double wavelength, double NA, double optical_power)
            {
                double omega = 0.61 * wavelength / NA;
                double area = 3.1415926535 * pow(omega / 2, 2);
                double intensity = optical_power / area;
                return sqrt(2.0 * intensity / (C_ * EPSILON0));
            }
        };
}
