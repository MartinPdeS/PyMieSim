#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.units import ureg


from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single import Setup



# Reference values from Scott Prahl's Mie scattering database
scott_prahl_values = {
    "Qsca": 1.1759,
    "Qext": 2.6257,
    "g": 0.80335,
}


@pytest.mark.parametrize(
    "measure_str", scott_prahl_values.keys(), ids=scott_prahl_values.keys()
)
def test_validation_scott_prahl(measure_str):
    """
    Validate PyMieSim results against known values from Scott Prahl's Mie scattering database.

    Parameters:
        measure_str (str): The scattering parameter to validate (e.g., 'Qsca', 'Qext', 'g').
    """
    polarization_state = PolarizationState(angle=0 * ureg.degree)

    source = Gaussian(
        wavelength=1000 * ureg.nanometer,
        polarization=polarization_state,
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )

    scatterer = Sphere(
        diameter=1000 * ureg.nanometer,
        material=1.5 + 0.5j * ureg.RIU,
        medium=1.0 * ureg.RIU,
    )

    setup = Setup(
        source=source,
        scatterer=scatterer,
    )

    # Retrieve the reference value from Scott Prahl's data
    scott_prahl_value = scott_prahl_values[measure_str]

    # Retrieve the corresponding PyMieSim value for the same scattering parameter
    pymiesim_value = setup.get(measure_str)

    # Compare the PyMieSim value to the reference value with a relative tolerance of 0.1%
    assert np.isclose(pymiesim_value, scott_prahl_value, atol=0, rtol=1e-3), (
        f"Mismatch in {measure_str}: PyMieSim value = {pymiesim_value}, "
        f"Scott Prahl value = {scott_prahl_value}"
    )


if __name__ == "__main__":
    pytest.main([__file__])
