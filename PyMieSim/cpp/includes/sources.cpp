#pragma once

#include "definitions.cpp"

namespace SOURCE
{
    struct State
    {
        double
            wavelength,
            amplitude,
            k;

        std::vector<complex128>
            jones_vector;

        size_t
            amplitude_idx;

        bool
            is_polarized = true;


        State(){}
        State(
            const double  &wavelength,
            const std::vector<complex128> &jones_vector,
            const double &amplitude)
        :
            wavelength(wavelength),
            jones_vector(jones_vector),
            amplitude(amplitude)
        {
            if (std::isnan(real(jones_vector[0])))
                is_polarized = false;

            this->k = 2.0 * PI / this->wavelength;
        }
    };

    class Source
    {
        public:
        double
            wavelength,
            polarization,
            amplitude,
            k;

        Source(
            const double &wavelength,
            const double &polarization,
            const double amplitude)
        :
            wavelength(wavelength),
            polarization(polarization),
            amplitude(amplitude)
        {
            k = 2 * PI / wavelength;
        }


        Source(){}
        ~Source(){}

        void set_wavelength(double value) { wavelength = value; k = 2.0 * PI / value; }

        void set_polarization(double value) { polarization = value; }

        void set_amplitude(double value) { amplitude = value; }
    };
}