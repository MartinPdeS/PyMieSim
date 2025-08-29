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
        double angular_frequency;
        std::vector<size_t> indices;
        size_t wavelength_index;

    BaseSource() = default;
    BaseSource(double wavelength, std::vector<complex128> jones_vector, double amplitude);

    protected:
        /**
         * @brief Computes the angular frequency from wavelength.
         *
         * Angular frequency ω = 2πc/λ, where c is the speed of light in vacuum.
         *
         * @param wavelength Wavelength in meters.
         * @return Angular frequency in rad/s.
         */
        double compute_angular_frequency(double wavelength) const;

        /**
         * @brief Computes the wavenumber from wavelength.
         *
         * Wavenumber k = 2π/λ.
         *
         * @param wavelength Wavelength in meters.
         * @return Wavenumber in rad/m.
         */
        double compute_wavenumber(double wavelength) const;

        /**
         * @brief Updates derived quantities (wavenumber, angular frequency) from wavelength.
         *
         * This method should be called whenever the wavelength is modified to ensure
         * all derived quantities remain consistent.
         */
        void update_derived_quantities();

};

class Planewave: public BaseSource {
    public:
        Planewave() = default;

        /**
         * @brief Constructs a plane wave source.
         *
         * @param wavelength Wavelength in meters.
         * @param jones_vector Polarization state as Jones vector [Ex, Ey].
         * @param amplitude Electric field amplitude.
         */
        Planewave(double wavelength, std::vector<complex128> jones_vector, double amplitude);
    };

class Gaussian: public BaseSource {
    public:
        double NA;
        double optical_power;

        Gaussian() = default;

        /**
         * @brief Constructs a Gaussian beam source.
         *
         * @param wavelength Wavelength in meters.
         * @param jones_vector Polarization state as Jones vector [Ex, Ey].
         * @param NA Numerical aperture of the Gaussian beam.
         * @param optical_power Optical power in watts.
         */
        Gaussian(double wavelength, std::vector<complex128> jones_vector, double NA, double optical_power);

        /**
         * @brief Computes the electric field amplitude from optical power.
         *
         * The amplitude is computed based on the relationship between optical power,
         * beam area, and electric field intensity for a Gaussian beam.
         *
         * @param wavelength Wavelength in meters.
         * @param NA Numerical aperture.
         * @param optical_power Optical power in watts.
         * @return Electric field amplitude.
         */
        double compute_amplitude_from_power(double wavelength, double NA, double optical_power);
    };
