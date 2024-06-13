#pragma once

#include <vector>
#include <complex>
#include <cmath> // For std::isnan

#define PI (double)3.14159265358979323846264338
#define EPSILON0 (double)8.854187817620389e-12
#define C_ (double)299792458.0
typedef std::complex<double> complex128;

namespace SOURCE
{
    class BaseSource {
        public:
            double wavelength = 0.0;
            std::vector<complex128> jones_vector;
            double amplitude = 0.0;
            double k = 0.0;

        BaseSource() = default;
        BaseSource(double wavelength, std::vector<complex128> jones_vector, double amplitude): wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude)
        {
            this->k = 2 * PI / this->wavelength;
        }

    };

    class Planewave: public BaseSource {
        public:
            Planewave() = default;

            Planewave(double wavelength, std::vector<complex128> jones_vector, double amplitude) {
                this->wavelength = wavelength;
                this->k = 2 * PI / this->wavelength;
                this->jones_vector = jones_vector;
                this->amplitude = amplitude;
            }
        };

    class Gaussian: public BaseSource {
        public:
            double NA = 0.0;
            double optical_power = 0.0;

            Gaussian() = default;

            Gaussian(double wavelength, std::vector<complex128> jones_vector, double NA, double optical_power)
            : NA(NA), optical_power(optical_power) {
                this->wavelength = wavelength;
                this->k = 2 * PI / this->wavelength;
                this->jones_vector = jones_vector;
                this->compute_amplitude_from_power();
            }

            void compute_amplitude_from_power()
            {
                double omega = 0.61 * wavelength / NA;
                double area = 3.1415926535 * pow(omega / 2, 2);
                double intensity = optical_power / area;
                this->amplitude = sqrt(2.0 * intensity / (C_ * EPSILON0));
            }
        };

    class Set
    {
        public:
            std::vector<double> wavelength;
            std::vector<std::vector<complex128>> jones_vector;
            std::vector<double> amplitude;
            std::vector<size_t> shape;

            Set() = default;

            Set(const std::vector<double> &wavelength,
                const std::vector<std::vector<complex128>> &jones_vector,
                const std::vector<double> &amplitude
            ) : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude)
            {
                this->shape = {this->wavelength.size(), this->jones_vector.size()};
            }

            Planewave to_object(size_t index_wavelength, size_t index_jones) const
            {
                return Planewave(
                    this->wavelength[index_wavelength],
                    this->jones_vector[index_jones],
                    this->amplitude[index_wavelength]
                );
            }
    };
}
