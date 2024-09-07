#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup, measure
import PyMieScatt as ps

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


def get_pymiesim_data(source, index, diameters, measure_string: str):
    """
    Retrieve simulation data using PyMieSim for a spherical scatterer.

    Parameters:
        source (Gaussian): The light source.
        index (complex): Refractive index of the sphere.
        diameters (array): Array of sphere diameters.
        measure_string (str): The type of measurement (e.g., 'Qext', 'Qsca').

    Returns:
        np.ndarray: Simulated data from PyMieSim.
    """
    scatterer = Sphere(
        diameter=diameters,
        index=index,
        medium_index=1.0,  # Assuming air as the medium
        source=source
    )

    experiment = Setup(
        scatterer=scatterer,
        source=source,
        detector=None
    )

    # Retrieve the specified measurement from the experiment
    data = experiment.get(getattr(measure, measure_string), export_as_numpy=True)
    return data.squeeze()


def get_pymiescatt_data(source, index, diameters, measure_string: str):
    """
    Retrieve simulation data using PyMieScatt for a spherical scatterer.

    Parameters:
        source (Gaussian): The light source.
        index (complex): Refractive index of the sphere.
        diameters (array): Array of sphere diameters.
        measure_string (str): The type of measurement (e.g., 'Qext', 'Qsca').

    Returns:
        np.ndarray: Simulated data from PyMieScatt.
    """
    pymiescatt_data = []
    for diameter in diameters:
        # Get scattering efficiencies from PyMieScatt
        efficiencies = ps.MieQ(
            m=index,
            wavelength=source.wavelength[0],
            diameter=diameter,
        )

        # Extract the required measurement
        measure_idx = PYMIESCATT_MEASURES.get(measure_string)
        pymiescatt_data.append(float(efficiencies[measure_idx]))

    return np.asarray(pymiescatt_data).squeeze()


def get_comparison(wavelength, index, diameters, measure_string: str):
    """
    Compare PyMieSim and PyMieScatt results for given parameters.

    Parameters:
        wavelength (float): Wavelength of the source in meters.
        index (complex): Refractive index of the sphere.
        diameters (array): Array of sphere diameters.
        measure_string (str): The type of measurement (e.g., 'Qext', 'Qsca').

    Returns:
        tuple: (PyMieSim results, PyMieScatt results)
    """
    # Create a Gaussian light source
    source = Gaussian(
        wavelength=wavelength,
        polarization=0,
        optical_power=1e-3,  # Optical power in Watts
        NA=0.2
    )

    # Get data from PyMieScatt
    pymiescatt_data = get_pymiescatt_data(
        source=source,
        index=index,
        diameters=diameters,
        measure_string=measure_string
    )

    # Get data from PyMieSim
    pymiesim_data = get_pymiesim_data(
        source=source,
        index=index,
        diameters=diameters,
        measure_string=measure_string
    )

    return pymiesim_data, pymiescatt_data


@pytest.mark.parametrize('measure_string', ['Qext', 'Qsca', 'Qabs', 'g', 'Qpr'])
def test_comparison(measure_string: str):
    """
    Test comparison between PyMieSim and PyMieScatt data for various scattering parameters.

    Parameters:
        measure_string (str): The type of measurement to compare (e.g., 'Qext', 'Qsca').
    """
    pymiesim_data, pymiescatt_data = get_comparison(
        wavelength=632e-9,  # Wavelength in meters (e.g., 632 nm)
        index=1.4 + 0.3j,   # Complex refractive index of the sphere
        diameters=np.geomspace(10e-9, 6000e-9, 800),  # Log-spaced array of diameters
        measure_string=measure_string
    )

    # Compare the data with a relative tolerance of 0.001 (0.1%)
    discrepancy = np.isclose(pymiesim_data, pymiescatt_data, atol=0, rtol=1e-3)

    # Ensure that more than 90% of the points are within the acceptable tolerance
    assert discrepancy.mean() > 0.9, f"Mismatch in PyMieSim vs PyMieScatt for {measure_string}"


if __name__ == "__main__":
    pytest.main([__file__])
