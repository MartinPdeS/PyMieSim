#include "source/source.h"


// ---------------------- Constructors ---------------------------------------

BaseSource::BaseSource(double wavelength, std::vector<complex128> jones_vector, double amplitude)
: wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude), wavenumber(2 * PI / wavelength)
{}

Planewave::Planewave(double wavelength, std::vector<complex128> jones_vector, double amplitude)
: BaseSource(wavelength, jones_vector, amplitude)
{}

Gaussian::Gaussian(double wavelength, std::vector<complex128> jones_vector, double NA, double optical_power)
: BaseSource(wavelength, jones_vector, compute_amplitude_from_power(wavelength, NA, optical_power)), NA(NA), optical_power(optical_power)
{}

// ---------------------- Others ---------------------------------------

double Gaussian::compute_amplitude_from_power(double wavelength, double NA, double optical_power)
{
    double omega = 0.61 * wavelength / NA;
    double area = 3.1415926535 * pow(omega / 2, 2);
    double intensity = optical_power / area;
    return sqrt(2.0 * intensity / (C_ * EPSILON0));
}

