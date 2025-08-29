#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import pandas as pd
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
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
    return Gaussian(
        wavelength=1000 * ureg.nanometer,  # Wavelength in meters
        polarization=0 * ureg.degree,  # Polarization angle
        optical_power=1 * ureg.watt,  # Optical power in ureg.watts
        NA=0.3 * ureg.AU,  # Numerical aperture
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
    core_index = 1.4 + 0.3j * ureg.RIU
    core_index = 1.4 + 0.3j * ureg.RIU
    shell_index = 1.3 * ureg.RIU
    core_diameters = np.geomspace(10, 6000, 80) * ureg.nanometer
    shell_thickness = 600 * ureg.nanometer
    medium_indexes = [1.0] * ureg.RIU

    # Get data from PyMieSim
    scatterer = CoreShell(
        core_diameter=core_diameters,
        core_property=core_index,
        shell_thickness=shell_thickness,
        shell_property=shell_index,
        medium_property=medium_indexes,
        source=gaussian_source,
    )

    experiment = Setup(scatterer=scatterer, source=gaussian_source)

    # Retrieve the specified measurement from the experiment
    pymiesim_data = experiment.get(measure, drop_unique_level=False)

    discrepency = np.allclose(
        pymiesim_data[measure].squeeze().values.quantity,
        pymiescatt_dataframe[measure].squeeze().values,
        atol=1e-6,
        rtol=1e-2,
    )

    assert discrepency, f"Mismatch in PyMieSim vs PyMieScatt for {measure}"


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
