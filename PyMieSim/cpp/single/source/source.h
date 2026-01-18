#pragma once

#include <array>
#include <complex>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "single/polarization/polarization.h"

#define PI (double)3.14159265358979323846264338
#define EPSILON0 (double)8.854187817620389e-12
#define C_ (double)299792458.0

using complex128 = std::complex<double>;

class BaseSource {
public:
    double wavelength = 0.0;

    // Replace raw Jones vector storage with a polarization object.
    JonesVector polarization;

    double amplitude = 0.0;
    double wavenumber = 0.0;
    double angular_frequency = 0.0;

    std::vector<std::size_t> indices;
    std::size_t wavelength_index = 0;

    BaseSource() = default;
    virtual ~BaseSource() = default;

    // New canonical constructor
    BaseSource(double wavelength_meter, JonesVector polarization_value, double amplitude_value);

    // Optional backward compatible constructor for existing call sites
    BaseSource(double wavelength_meter, std::vector<complex128> jones_vector, double amplitude_value);

    // Convenience: first row Jones vector as [Ex, Ey]
    std::array<complex128, 2> get_jones_vector_first_row() const
    {
        const auto& rows = this->polarization.elements();
        if (rows.empty())
            throw std::runtime_error("polarization has no elements.");

        return rows[0];
    }

protected:
    double compute_angular_frequency(double wavelength_meter) const;

    double compute_wavenumber(double wavelength_meter) const;

    void update_derived_quantities();

    void validate_polarization() const
    {
        // In your Python model, element is (N, 2). We enforce at least one row.
        const auto& rows = this->polarization.elements();

        if (rows.empty()) {
            throw std::invalid_argument("polarization must contain at least one Jones vector row.");
        }

        // Also enforce finite values
        const auto v = rows[0];
        if (!std::isfinite(v[0].real()) || !std::isfinite(v[0].imag()) ||
            !std::isfinite(v[1].real()) || !std::isfinite(v[1].imag())) {
            throw std::invalid_argument("polarization contains non finite values.");
        }
    }

    static JonesVector jones_vector_to_polarization(std::vector<complex128> jones_vector)
    {
        if (jones_vector.size() != 2) {
            throw std::invalid_argument("jones_vector must have exactly 2 components [Ex, Ey].");
        }

        JonesVector::Row row{jones_vector[0], jones_vector[1]};
        return JonesVector(row);
    }
};

class Planewave : public BaseSource {
public:
    Planewave() = default;

    // New canonical constructor
    Planewave(
        double wavelength_meter,
        JonesVector polarization_value,
        double amplitude_value
    );

    // Optional backward compatible constructor
    Planewave(
        double wavelength_meter,
        std::vector<complex128> jones_vector,
        double amplitude_value
    );
};

class Gaussian : public BaseSource {
public:
    double NA = 0.0;
    double optical_power = 0.0;
    double waist = 0.0;
    double peak_intensity = 0.0;
    double area = 0.0;

    Gaussian() = default;

    Gaussian(
        double wavelength_meter,
        JonesVector polarization_value,
        double NA_value,
        double optical_power_watt
    );

    Gaussian(
        double wavelength_meter,
        std::vector<complex128> jones_vector,
        double NA_value,
        double optical_power_watt
    );

    /*
    @brief Compute the amplitude from optical power, wavelength, and NA.
    @param wavelength_meter Wavelength in meters.
    @param NA_value Numerical aperture (dimensionless).
    @param optical_power_watt Optical power in watts.
    @return Amplitude in V/m.
    @throws std::invalid_argument if any parameter is non-positive.
    */
    double compute_amplitude_from_power(double wavelength_meter, double NA_value, double optical_power_watt);
};
