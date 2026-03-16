#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import pandas as pd
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import CoreShellSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.directories import validation_data_path


# Mapping of measurement types to PyMieScatt indices
PYMIESCATT_MEASURES = {
    "Qext": 0,
    "Qsca": 1,
    "Qabs": 2,
    "g": 3,
    "Qpr": 4,
    "Qback": 5,
    "Qratio": 6,
}


@pytest.fixture
def gaussian_source():
    return GaussianSet(
        wavelength=1000 * ureg.nanometer,
        polarization=PolarizationSet(angles=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3,
    )


@pytest.fixture
def pymiescatt_dataframe():
    filename = validation_data_path / "pymiescatt/validation_coreshell.csv"
    return pd.read_csv(filename)


@pytest.mark.parametrize("measure", ["Qext", "Qsca", "Qabs", "g", "Qpr"])
def test_comparison(pymiescatt_dataframe, gaussian_source, measure: str):
    """
    Test comparison between PyMieSim and PyMieScatt data for various scattering parameters.

    Parameters
    ----------
    measure : str
        The type of measurement to compare (e.g., 'Qext', 'Qsca').
    """
    core_index = [1.4 + 0.3j]
    core_index = [1.4 + 0.3j]
    shell_index = [1.3]
    core_diameters = np.geomspace(10, 6000, 80) * ureg.nanometer
    shell_thickness = [300] * ureg.nanometer
    medium_indexes = [1.0]

    # Get data from PyMieSim
    scatterer = CoreShellSet(
        core_diameter=core_diameters,
        core_material=core_index,
        shell_thickness=shell_thickness,
        shell_material=shell_index,
        medium=medium_indexes,
    )

    experiment = Setup(scatterer_set=scatterer, source_set=gaussian_source)

    # Retrieve the specified measurement from the experiment
    pymiesim_data = experiment.get(measure, drop_unique_level=True)

    discrepency = np.allclose(
        pymiesim_data.get(measure).magnitude,
        pymiescatt_dataframe[measure].squeeze().values,
        atol=1e-6,
        rtol=1e-2,
    )

    assert discrepency, f"Mismatch in PyMieSim vs PyMieScatt for {measure}"


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
