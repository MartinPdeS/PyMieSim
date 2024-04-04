#pragma once

#include "definitions.cpp" // Ensure this file provides the necessary definitions like PI and complex128
#include <vector>
#include <complex>
#include <cmath> // For std::isnan

namespace SOURCE {
    struct State {
        double wavelength = 0.0;
        double amplitude = 0.0;
        double k = 0.0;
        std::vector<std::complex<double>> jones_vector;
        size_t amplitude_idx = 0;
        bool is_polarized = true;

        // Default constructor
        State() = default;

        // Parameterized constructor
        State(double wavelength, const std::vector<std::complex<double>>& jones_vector, double amplitude)
            : wavelength(wavelength), amplitude(amplitude), jones_vector(jones_vector) {
            if (jones_vector.empty() || std::isnan(std::real(jones_vector[0]))) {
                is_polarized = false;
            }
            this->k = 2.0 * PI / this->wavelength;
        }
    };

    class Source {
    public:
        double wavelength = 0.0;
        double polarization = 0.0;
        double amplitude = 0.0;
        double k = 0.0;

        // Default constructor
        Source() = default;

        // Parameterized constructor
        Source(double wavelength, double polarization, double amplitude)
            : wavelength(wavelength), polarization(polarization), amplitude(amplitude) {
            this->k = 2 * PI / this->wavelength;
        }

        void set_wavelength(double value) {
            wavelength = value;
            k = 2.0 * PI / value;
        }

        void set_polarization(double value) {
            polarization = value;
        }

        void set_amplitude(double value) {
            amplitude = value;
        }
    };
}
