#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import pandas as pd
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.directories import validation_data_path


@pytest.fixture
def gaussian_source():
    return GaussianSet(
        wavelength=[1000] * ureg.nanometer,
        polarization=PolarizationSet(angles=0 * ureg.degree),
        optical_power=[1] * ureg.watt,
        numerical_aperture=[0.3] * ureg.AU,
    )


@pytest.fixture
def pymiescatt_dataframe():
    filename = validation_data_path / "pymiescatt/validation_sphere.csv"
    return pd.read_csv(filename)


@pytest.mark.parametrize("measure", ["Qext", "Qsca", "Qabs", "g", "Qpr"])
def test_comparison(pymiescatt_dataframe, gaussian_source, measure: str):
    """
    Test comparison between PyMieSim and PyMieScatt data for various scattering parameters.

    Parameters:
        measure (str): The type of measurement to compare (e.g., 'Qext', 'Qsca').
    """
    diameters = np.geomspace(10, 6000, 800) * ureg.nanometer

    # Get data from PyMieSim
    scatterer = SphereSet(
        diameter=diameters,
        material=[1.4 + 0.3j] * ureg.RIU,
        medium=[1.0] * ureg.RIU,
    )

    experiment = Setup(scatterer_set=scatterer, source_set=gaussian_source)

    # Retrieve the specified measurement from the experiment
    pymiesim_data = experiment.get(measure)

    discrepency = np.allclose(
        pymiesim_data.get(measure).magnitude,
        pymiescatt_dataframe[measure].squeeze().values,
        atol=1e-3,
        rtol=1e-2,
    )

    assert (
        discrepency
    ), f"Mismatch in PyMieSim vs PyMieScatt for {measure} \n {pymiesim_data[measure].squeeze().values.quantity}"


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
