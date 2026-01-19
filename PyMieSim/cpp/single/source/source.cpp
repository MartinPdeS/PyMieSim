#include "./source.h"


// ---------------------- BaseSource Implementation ---------------------------------------

BaseSource::BaseSource(
    double _wavelength,
    JonesVector polarization_value,
    double _amplitude
)
    : wavelength(_wavelength),
      polarization(std::move(polarization_value)),
      amplitude(_amplitude)
{
    this->update_derived_quantities();
    this->validate_polarization();
}

BaseSource::BaseSource(
    double _wavelength,
    std::vector<complex128> _jones_vector,
    double _amplitude
)
    : BaseSource(
        _wavelength,
        jones_vector_to_polarization(std::move(_jones_vector)),
        _amplitude
    )
{
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
    JonesVector polarization_value,
    double _amplitude)
    : BaseSource(_wavelength, std::move(polarization_value), _amplitude)
{
}

Planewave::Planewave(
    double _wavelength,
    std::vector<complex128> _jones_vector,
    double _amplitude)
    : BaseSource(_wavelength, std::move(_jones_vector), _amplitude)
{
}

// ---------------------- Gaussian Implementation ---------------------------------------

Gaussian::Gaussian(
    double _wavelength,
    JonesVector polarization_value,
    double _NA,
    double _optical_power)
    : BaseSource(_wavelength, std::move(polarization_value), 0.0),
      NA(_NA),
      optical_power(_optical_power)
{
    this->amplitude = this->compute_amplitude_from_power(this->wavelength, this->NA, this->optical_power);
}

Gaussian::Gaussian(
    double _wavelength,
    std::vector<complex128> _jones_vector,
    double _NA,
    double _optical_power)
    : Gaussian(
        _wavelength,
        jones_vector_to_polarization(std::move(_jones_vector)),
        _NA,
        _optical_power
    )
{
}

double Gaussian::compute_amplitude_from_power(
    double wavelength_meter,
    double NA_value,
    double optical_power_watt)
{
    if (wavelength_meter <= 0.0) {
        throw std::invalid_argument("wavelength must be positive.");
    }
    if (NA_value <= 0.0) {
        throw std::invalid_argument("NA must be positive.");
    }
    if (optical_power_watt <= 0.0) {
        throw std::invalid_argument("optical_power must be positive.");
    }

    this->waist = 0.61 * wavelength_meter / NA_value;
    this->peak_intensity = 2.0 * optical_power_watt / (Constants::PI * this->waist * this->waist);
    return std::sqrt( 2.0 * this->peak_intensity / (Constants::LIGHT_SPEED * Constants::EPSILON0) );
}
