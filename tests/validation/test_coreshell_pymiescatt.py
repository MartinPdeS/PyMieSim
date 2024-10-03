#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
import PyMieScatt as ps
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyMieSim.utils import get_pymiescatt_coreshell_dataframe

# Mapping of measurement types to PyMieScatt indices
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


@pytest.mark.parametrize('measure', ['Qext', 'Qsca', 'Qabs', 'g', 'Qpr'])
def test_comparison(gaussian_source, measure: str):
    """
    Test comparison between PyMieSim and PyMieScatt data for various scattering parameters.

    Parameters:
        measure (str): The type of measurement to compare (e.g., 'Qext', 'Qsca').
    """
    wavelength = gaussian_source.wavelength
    core_index = 1.4 + 0.3j * RIU
    core_index = 1.4 + 0.3j * RIU
    shell_index = 1.3 * RIU
    core_diameters = np.geomspace(10, 6000, 80) * nanometer
    shell_width = 600 * nanometer
    medium_indexes = [1.0] * RIU

    pymiescatt_df = get_pymiescatt_coreshell_dataframe(
        wavelengths=wavelength,
        shell_indexes=shell_index,
        core_indexes=core_index,
        core_diameters=core_diameters,
        shell_widths=shell_width,
        medium_indexes=medium_indexes
    )

    # Get data from PyMieSim
    scatterer = CoreShell(
        core_diameter=core_diameters,
        core_property=core_index,
        shell_width=shell_width,
        shell_property=shell_index,
        medium_property=medium_indexes,
        source=gaussian_source
    )

    experiment = Setup(scatterer=scatterer, source=gaussian_source)

    # Retrieve the specified measurement from the experiment
    pymiesim_data = experiment.get(measure, drop_unique_level=False)

    discrepency = np.allclose(
        pymiesim_data[measure].squeeze().values.quantity,
        pymiescatt_df[measure].squeeze().values.quantity,
        atol=1e-6,
        rtol=1e-2
    )

    assert discrepency, f"Mismatch in PyMieSim vs PyMieScatt for {measure}"


if __name__ == "__main__":
    pytest.main([__file__])
