#include "source/source.h"


// ---------------------- BaseSource Implementation ---------------------------------------

BaseSource::BaseSource(double wavelength, std::vector<complex128> jones_vector, double amplitude)
: wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude)
{
    update_derived_quantities();
}

double BaseSource::compute_angular_frequency(double wavelength) const {
    return 2.0 * PI * C_ / wavelength;
}

double BaseSource::compute_wavenumber(double wavelength) const {
    return 2.0 * PI / wavelength;
}

void BaseSource::update_derived_quantities() {
    wavenumber = compute_wavenumber(wavelength);
    angular_frequency = compute_angular_frequency(wavelength);
}

// ---------------------- Planewave Implementation ---------------------------------------

Planewave::Planewave(double wavelength, std::vector<complex128> jones_vector, double amplitude)
: BaseSource(wavelength, jones_vector, amplitude)
{}

// ---------------------- Gaussian Implementation ---------------------------------------

Gaussian::Gaussian(double wavelength, std::vector<complex128> jones_vector, double NA, double optical_power)
: BaseSource(wavelength, jones_vector, compute_amplitude_from_power(wavelength, NA, optical_power)),
  NA(NA), optical_power(optical_power)
{}

double Gaussian::compute_amplitude_from_power(double wavelength, double NA, double optical_power)
{
    double omega = 0.61 * wavelength / NA;
    double area = PI * pow(omega / 2, 2);
    double intensity = optical_power / area;
    return sqrt(2.0 * intensity / (C_ * EPSILON0));
}

