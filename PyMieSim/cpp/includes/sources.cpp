#pragma once

#include "definitions.cpp" // Ensure this file provides the necessary definitions like PI and complex128
#include <vector>
#include <complex>
#include <cmath> // For std::isnan
#include "sources.h"

namespace SOURCE {

    void Source::set_wavelength(double value) {
        wavelength = value;
        k = 2.0 * PI / value;
    }

    void Source::set_polarization(double value) {
        polarization = value;
    }

    void Source::set_amplitude(double value) {
        amplitude = value;
    }
}
