#pragma once

#include <array>
#include <complex>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "polarization/polarization.h"
#include <utils/constants.h>

using complex128 = std::complex<double>;

class BaseSource {
public:
    double wavelength = 0.0;

    // Replace raw Jones vector storage with a polarization object.
    PolarizationState polarization;

    double amplitude = 0.0;
    double wavenumber_vacuum = 0.0;
    double angular_frequency = 0.0;

    std::vector<std::size_t> indices;
    std::size_t wavelength_index = 0;

    BaseSource() = default;
    virtual ~BaseSource() = default;

    // New canonical constructor
    BaseSource(double wavelength_meter, PolarizationState polarization_value, double amplitude_value);

    // Optional backward compatible constructor for existing call sites
    BaseSource(double wavelength_meter, std::vector<complex128> jones_vector, double amplitude_value);

    // Convenience: first row Jones vector as [Ex, Ey]
    std::vector<complex128> get_jones_vector_first_row() const
    {
        return this->polarization.jones_vector;
    }

protected:
    double compute_angular_frequency(double wavelength_meter) const;

    double compute_wavenumber_vacuum(double wavelength_meter) const;

    void update_derived_quantities();

    void validate_polarization() const
    {
        auto row = this->polarization.jones_vector;
        if (!std::isfinite(row[0].real()) || !std::isfinite(row[0].imag()) ||
            !std::isfinite(row[1].real()) || !std::isfinite(row[1].imag())) {
            throw std::invalid_argument("polarization contains non finite values.");
        }
    }
};

class Planewave : public BaseSource {
public:
    Planewave() = default;

    // New canonical constructor
    Planewave(
        double wavelength_meter,
        PolarizationState polarization_state,
        double amplitude_value
    );

};

class Gaussian : public BaseSource {
public:
    double numerical_aperture = 0.0;
    double optical_power = 0.0;
    double waist = 0.0;
    double peak_intensity = 0.0;
    double area = 0.0;

    Gaussian() = default;

    Gaussian(
        double wavelength_meter,
        PolarizationState polarization_state,
        double numerical_aperture_value,
        double optical_power_watt
    );

    /*
    @brief Compute the amplitude from optical power, wavelength, and numerical aperture.
    @param wavelength_meter Wavelength in meters.
    @param numerical_aperture_value Numerical aperture (dimensionless).
    @param optical_power_watt Optical power in watts.
    @return Amplitude in V/m.
    @throws std::invalid_argument if any parameter is non-positive.
    */
    double compute_amplitude_from_power(double wavelength_meter, double numerical_aperture_value, double optical_power_watt);
};
