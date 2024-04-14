#pragma once

#include "definitions.cpp" // Ensure this file provides the necessary definitions like PI and complex128
#include <vector>
#include <complex>
#include <cmath> // For std::isnan

namespace SOURCE
{
    class Set
    {
        public:
            std::vector<double> wavelength;
            std::vector<double> amplitude;
            std::vector<std::vector<complex128>> jones_vector;
            std::vector<size_t> shape;

            Set() = default;

            Set(
                const std::vector<double> &wavelength,
                const std::vector<std::vector<complex128>> &jones_vector,
                const std::vector<double> &amplitude
            ) : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude)
            {
                this->shape = {this->wavelength.size(), this->jones_vector.size()};
            }
    };

    class BaseSource {
        public:
            double wavelength = 0.0;
            std::vector<std::complex<double>> jones_vector;
            double amplitude = 0.0;
            double k = 0.0;

        BaseSource() = default;
        BaseSource(double wavelength, std::vector<std::complex<double>> jones_vector, double amplitude): wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude)
        {
            this->k = 2 * PI / this->wavelength;
        }

    };

    class Planewave: public BaseSource {
        public:
            Planewave() = default;

            Planewave(double wavelength, std::vector<std::complex<double>> jones_vector, double amplitude)
            : BaseSource(wavelength, jones_vector, amplitude) {}
        };

    class Gaussian: public BaseSource {
        public:
            double NA = 0.0;
            double optical_power = 0.0;

            Gaussian() = default;

            Gaussian(double wavelength, std::vector<std::complex<double>> jones_vector, double NA, double optical_power) 
            : NA(NA), optical_power(optical_power) {
                this->compute_amplitude_from_power();
                BaseSource(wavelength, jones_vector, amplitude);

            }

            double compute_amplitude_from_power()
            {
                double omega = 0.61 * wavelength / NA;
                double area = 3.1415926535 * pow(omega / 2, 2);
                double intensity = optical_power / area;
                double epsilon_0 = 8.8541878128e-12;
                double c = 299792458.0;
                this->amplitude = sqrt(2.0 * intensity / (c * epsilon_0));
            }
        };
}
