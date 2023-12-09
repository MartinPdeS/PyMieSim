#pragma once

#include "definitions.cpp"

namespace SOURCE
{

  struct State
  {
    double wavelength;
    CVector jones_vector;
    double amplitude;
    size_t amplitude_idx;
    double k;
    bool is_polarized = true;


    State(){}
    State(double  &wavelength, CVector &jones_vector, double &amplitude)
    : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude)
    {
      if (std::isnan( real( jones_vector[0] ))) is_polarized = false;

      this->k = 2.0 * PI / this->wavelength;
    }
  };

  class Source
  {
  public:
    double wavelength, polarization, amplitude, k;

    Source(double &wavelength, double &polarization, double amplitude)
    : wavelength(wavelength), polarization(polarization), amplitude(amplitude)
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