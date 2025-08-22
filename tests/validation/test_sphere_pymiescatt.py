#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import pandas as pd
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.directories import validation_data_path


@pytest.fixture
def gaussian_source():
    return Gaussian(
        wavelength=1000 * ureg.nanometer,  # Wavelength in meters
        polarization=0 * ureg.degree,   # Polarization angle
        optical_power=1 * ureg.watt,  # Optical power in ureg.watts
        NA=0.3 * ureg.AU           # Numerical aperture
    )


@pytest.fixture
def pymiescatt_dataframe():
    filename = validation_data_path / "pymiescatt/validation_sphere.csv"
    return pd.read_csv(filename)


@pytest.mark.parametrize('measure',  ['Qext', 'Qsca', 'Qabs', 'g', 'Qpr'])
def test_comparison(pymiescatt_dataframe, gaussian_source, measure: str):
    """
    Test comparison between PyMieSim and PyMieScatt data for various scattering parameters.

    Parameters:
        measure (str): The type of measurement to compare (e.g., 'Qext', 'Qsca').
    """
    index = (1.4 + 0.3j) * ureg.RIU
    diameters = np.geomspace(10, 6000, 800) * ureg.nanometer
    medium_indexes = [1.0] * ureg.RIU

    # Get data from PyMieSim
    scatterer = Sphere(
        diameter=diameters,
        property=index,
        medium_property=medium_indexes,
        source=gaussian_source
    )

    experiment = Setup(scatterer=scatterer, source=gaussian_source)

    # Retrieve the specified measurement from the experiment
    pymiesim_data = experiment.get(measure, drop_unique_level=False)

    discrepency = np.allclose(
        pymiesim_data[measure].squeeze().values.quantity,
        pymiescatt_dataframe[measure].squeeze().values,
        atol=1e-3,
        rtol=1e-2
    )

    assert discrepency, f"Mismatch in PyMieSim vs PyMieScatt for {measure} \n {pymiesim_data[measure].squeeze().values.quantity}"


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
