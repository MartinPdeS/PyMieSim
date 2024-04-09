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

        State() = default;

        State(double wavelength, const std::vector<std::complex<double>>& jones_vector, double amplitude)
            : wavelength(wavelength), amplitude(amplitude), jones_vector(jones_vector) {
            if (jones_vector.empty() || std::isnan(std::real(jones_vector[0]))) {
                is_polarized = false;
            }
            this->k = 2.0 * PI / this->wavelength;
        }
    };

    class Set
    {
        public:
            std::vector<double> wavelength;
            std::vector<double> amplitude;
            std::vector<std::vector<complex128>> jones_vector;
            std::vector<State> state_list;

            Set() = default;

            Set(const std::vector<double> &wavelength,
                const std::vector<std::vector<complex128>> &jones_vector,
                const std::vector<double> &amplitude)
            : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude){}

            State operator[](const size_t &idx){return this->state_list[idx];}

            std::vector<size_t> get_array_shape() const
            {
                return {this->wavelength.size(), this->jones_vector.size()};
            }

            size_t get_array_size() const
            {
                std::vector<size_t> full_shape = this->get_array_shape();

                size_t full_size = 1;

                for (auto e: full_shape)
                    full_size *= e;

                return full_size;
            }
    };

    class Source {
    public:
        double wavelength = 0.0;
        double polarization = 0.0;
        double amplitude = 0.0;
        double k = 0.0;

        Source() = default;

        Source(double wavelength, double polarization, double amplitude)
            : wavelength(wavelength), polarization(polarization), amplitude(amplitude) {
            this->k = 2 * PI / this->wavelength;
        }

        void set_wavelength(double value);
        void set_polarization(double value);
        void set_amplitude(double value);
    };
}
