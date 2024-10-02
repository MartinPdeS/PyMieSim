#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup, measure
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyMieSim.utils import get_pymiescatt_sphere_dataframe

# Mapping PyMieScatt measure strings to their respective indices
PYMIESCATT_MEASURES = {
    'Qext': 0,
    'Qsca': 1,
    'Qabs': 2,
    'g': 3,
    'Qpr': 4,
    'Qback': 5,
    'Qratio': 6
}

@pytest.fixture
def gaussian_source():
    return Gaussian(
        wavelength=1000 * nanometer,  # Wavelength in meters
        polarization=0 * degree,   # Polarization angle
        optical_power=1 * watt,  # Optical power in watts
        NA=0.3 * AU           # Numerical aperture
    )


@pytest.mark.parametrize('measure_string', ['Qext', 'Qsca', 'Qabs', 'g', 'Qpr'])
def test_comparison(gaussian_source, measure_string: str):
    """
    Test comparison between PyMieSim and PyMieScatt data for various scattering parameters.

    Parameters:
        measure_string (str): The type of measurement to compare (e.g., 'Qext', 'Qsca').
    """
    wavelength = 632 * nanometer
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
        diameter=diameters * 1.58231,
        index=index,
        medium_index=medium_indexes,
        source=gaussian_source
    )

    experiment = Setup(scatterer=scatterer, source=gaussian_source)

    # Retrieve the specified measurement from the experiment
    pymiesim_data = experiment.get(measure_string)

    discrepency = np.allclose(
        pymiesim_data[measure_string].values.numpy_data,
        pymiescatt_df[measure_string].values.numpy_data,
        atol=1e-6,
        rtol=1e-2
    )

    assert discrepency, f"Mismatch in PyMieSim vs PyMieScatt for {measure_string}"


if __name__ == "__main__":
    pytest.main([__file__])
