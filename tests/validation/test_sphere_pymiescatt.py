#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyMieSim.utils import get_pymiescatt_sphere_dataframe


@pytest.fixture
def gaussian_source():
    return Gaussian(
        wavelength=1000 * nanometer,  # Wavelength in meters
        polarization=0 * degree,   # Polarization angle
        optical_power=1 * watt,  # Optical power in watts
        NA=0.3 * AU           # Numerical aperture
    )


@pytest.mark.parametrize('measure',  ['Qext', 'Qsca', 'Qabs', 'g', 'Qpr'])
def test_comparison(gaussian_source, measure: str):
    """
    Test comparison between PyMieSim and PyMieScatt data for various scattering parameters.

    Parameters:
        measure (str): The type of measurement to compare (e.g., 'Qext', 'Qsca').
    """
    wavelength = gaussian_source.wavelength
    index = (1.4 + 0.3j) * RIU
    diameters = np.geomspace(10, 6000, 800) * nanometer
    medium_indexes = [1.0] * RIU

    pymiescatt_df = get_pymiescatt_sphere_dataframe(
        wavelengths=wavelength,
        indexes=index,
        diameters=diameters,
        medium_indexes=medium_indexes
    )

    # Get data from PyMieSim
    scatterer = Sphere(
        diameter=diameters,
        index=index,
        medium_index=medium_indexes,
        source=gaussian_source
    )

    experiment = Setup(scatterer=scatterer, source=gaussian_source)

    # Retrieve the specified measurement from the experiment
    pymiesim_data = experiment.get(measure, drop_unique_level=False)

    discrepency = np.allclose(
        pymiesim_data[measure].squeeze().values.quantity,
        pymiescatt_df[measure].squeeze().values.quantity,
        atol=1e-3,
        rtol=1e-2
    )

    assert discrepency, f"Mismatch in PyMieSim vs PyMieScatt for {measure}"


if __name__ == "__main__":
    pytest.main([__file__])
