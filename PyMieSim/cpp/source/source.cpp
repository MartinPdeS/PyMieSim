#include "./source.h"


// ---------------------- BaseSource Implementation ---------------------------------------

BaseSource::BaseSource(double _wavelength, std::vector<complex128> _jones_vector, double _amplitude)
: wavelength(_wavelength), jones_vector(_jones_vector), amplitude(_amplitude)
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

Planewave::Planewave(double _wavelength, std::vector<complex128> _jones_vector, double _amplitude)
: BaseSource(_wavelength, _jones_vector, _amplitude)
{}

// ---------------------- Gaussian Implementation ---------------------------------------

Gaussian::Gaussian(double _wavelength, std::vector<complex128> _jones_vector, double _NA, double _optical_power)
: BaseSource(_wavelength, _jones_vector, compute_amplitude_from_power(_wavelength, _NA, _optical_power)),
  NA(_NA), optical_power(_optical_power)
{}

double Gaussian::compute_amplitude_from_power(double wavelength, double NA, double optical_power)
{
    double omega = 0.61 * wavelength / NA;
    double area = PI * pow(omega / 2, 2);
    double intensity = optical_power / area;
    return sqrt(2.0 * intensity / (C_ * EPSILON0));
}
