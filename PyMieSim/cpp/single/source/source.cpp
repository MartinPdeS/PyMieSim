#include "./source.h"


// ---------------------- BaseSource Implementation ---------------------------------------

BaseSource::BaseSource(
    double _wavelength,
    PolarizationState _polarization_state,
    double _amplitude
)
    : wavelength(_wavelength),
      polarization(_polarization_state),
      amplitude(_amplitude)
{
    this->update_derived_quantities();
    this->validate_polarization();
}


double BaseSource::compute_angular_frequency(double wavelength_meter) const
{
    if (wavelength_meter <= 0.0) {
        throw std::invalid_argument("wavelength must be positive.");
    }
    return 2.0 * Constants::PI * Constants::LIGHT_SPEED / wavelength_meter;
}

double BaseSource::compute_wavenumber_vacuum(double wavelength_meter) const
{
    if (wavelength_meter <= 0.0) {
        throw std::invalid_argument("wavelength must be positive.");
    }
    return 2.0 * Constants::PI / wavelength_meter;
}

void BaseSource::update_derived_quantities()
{
    this->wavenumber_vacuum = this->compute_wavenumber_vacuum(this->wavelength);
    this->angular_frequency = this->compute_angular_frequency(this->wavelength);
}

// ---------------------- Planewave Implementation ---------------------------------------

Planewave::Planewave(
    double _wavelength,
    PolarizationState polarization_state,
    double _amplitude)
    : BaseSource(_wavelength, polarization_state, _amplitude)
{
}


// ---------------------- Gaussian Implementation ---------------------------------------

Gaussian::Gaussian(
    double _wavelength,
    PolarizationState polarization_state,
    double _numerical_aperture,
    double _optical_power)
    : BaseSource(_wavelength, polarization_state, 0.0),
      numerical_aperture(_numerical_aperture),
      optical_power(_optical_power)
{
    this->amplitude = this->compute_amplitude_from_power(this->wavelength, this->numerical_aperture, this->optical_power);
}


double Gaussian::compute_amplitude_from_power(
    double wavelength_meter,
    double numerical_aperture_value,
    double optical_power_watt)
{
    if (wavelength_meter <= 0.0) {
        throw std::invalid_argument("wavelength must be positive.");
    }
    if (numerical_aperture_value <= 0.0) {
        throw std::invalid_argument("Numerical aperture must be positive.");
    }
    if (optical_power_watt <= 0.0) {
        throw std::invalid_argument("optical_power must be positive.");
    }

    this->waist = 0.61 * wavelength_meter / numerical_aperture_value;
    this->peak_intensity = 2.0 * optical_power_watt / (Constants::PI * this->waist * this->waist);
    return std::sqrt( 2.0 * this->peak_intensity / (Constants::LIGHT_SPEED * Constants::EPSILON0) );
}
