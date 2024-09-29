#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# Reference values from Scott Prahl's Mie scattering database
scott_prahl_values = {
    'Qsca': 1.1759,  # Scattering efficiency
    'Qext': 2.6257,  # Extinction efficiency
    'g': 0.80335,    # Asymmetry parameter
}

@pytest.mark.parametrize('measure_str', scott_prahl_values.keys(), ids=scott_prahl_values.keys())
def test_validation_scott_prahl(measure_str):
    """
    Validate PyMieSim results against known values from Scott Prahl's Mie scattering database.

    Parameters:
        measure_str (str): The scattering property to validate (e.g., 'Qsca', 'Qext', 'g').
    """
    # Create a Gaussian light source
    source = Gaussian(
        wavelength=1000 * nanometer,   # Wavelength in meters (e.g., 1 micron)
        polarization=0 * degree,    # Polarization angle
        optical_power=1 * watt,   # Optical power in watts
        NA=0.3 * AU             # Numerical aperture
    )

    # Create a spherical scatterer
    scatterer = Sphere(
        diameter=1000 * nanometer,       # Diameter in meters (e.g., 1 micron)
        index=1.5 + 0.5j * RIU,    # Complex refractive index
        source=source,       # Associated light source
        medium_index=1.0 * RIU     # Refractive index of the medium (e.g., air)
    )

    # Retrieve the reference value from Scott Prahl's data
    scott_prahl_value = scott_prahl_values[measure_str]
    # Retrieve the corresponding PyMieSim value for the same scattering property
    pymiesim_value = getattr(scatterer, measure_str)

    # Compare the PyMieSim value to the reference value with a relative tolerance of 0.1%
    assert np.isclose(pymiesim_value, scott_prahl_value, atol=0, rtol=1e-3), (
        f'Mismatch in {measure_str}: PyMieSim value = {pymiesim_value}, '
        f'Scott Prahl value = {scott_prahl_value}'
    )


if __name__ == "__main__":
    pytest.main([__file__])
