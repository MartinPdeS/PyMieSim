#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure
import PyMieScatt as ps
from PyMieSim.units import nanometer, degree, watt, AU, RIU

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

def get_pymiesim_data(source, core_index, shell_index, core_diameters, shell_width, measure_string):
    """
    Retrieve simulation data using PyMieSim for a core-shell scatterer setup.

    Parameters:
        source (Gaussian): The light source.
        core_index (complex): Refractive index of the core.
        shell_index (float): Refractive index of the shell.
        core_diameters (array): Array of core diameters.
        shell_width (float): Shell thickness.
        measure_string (str): The type of measurement (e.g., 'Qext', 'Qsca').

    Returns:
        numpy.ndarray: Simulation results as a squeezed numpy array.
    """
    scatterer = CoreShell(
        core_diameter=core_diameters,
        shell_width=shell_width,
        core_index=core_index,
        shell_index=shell_index,
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


def get_pymiescatt_data(source, core_index, shell_index, core_diameters, shell_width, measure_string):
    """
    Retrieve simulation data using PyMieScatt for a core-shell scatterer setup.

    Parameters:
        source (Gaussian): The light source.
        core_index (complex): Refractive index of the core.
        shell_index (float): Refractive index of the shell.
        core_diameters (array): Array of core diameters.
        shell_width (float): Shell thickness.
        measure_string (str): The type of measurement (e.g., 'Qext', 'Qsca').

    Returns:
        numpy.ndarray: Simulation results as a squeezed numpy array.
    """
    pymiescatt_data = []
    for core_diameter in core_diameters:
        # Get Mie scattering efficiencies from PyMieScatt
        efficiencies = ps.MieQCoreShell(
            mCore=core_index,
            mShell=shell_index,
            wavelength=source.wavelength[0],
            dCore=core_diameter,
            dShell=core_diameter + shell_width
        )

        # Extract the required measurement
        measure_idx = PYMIESCATT_MEASURES.get(measure_string)
        pymiescatt_data.append(float(efficiencies[measure_idx]))

    return np.asarray(pymiescatt_data).squeeze()


def get_comparison(wavelength, core_index, shell_index, core_diameters, shell_width, measure_string):
    """
    Compare PyMieSim and PyMieScatt results for given parameters.

    Parameters:
        wavelength (float): Wavelength of the source.
        core_index (complex): Refractive index of the core.
        shell_index (float): Refractive index of the shell.
        core_diameters (array): Array of core diameters.
        shell_width (float): Shell thickness.
        measure_string (str): The type of measurement (e.g., 'Qext', 'Qsca').

    Returns:
        tuple: (PyMieSim results, PyMieScatt results)
    """
    # Create a Gaussian light source
    source = Gaussian(
        wavelength=wavelength,
        polarization=0 * degree,
        optical_power=1e-3 * watt,  # In Watts
        NA=0.2 * AU
    )

    # Get data from PyMieScatt
    pymiescatt_data = get_pymiescatt_data(
        source=source,
        core_index=core_index.magnitude,
        shell_index=shell_index.magnitude,
        core_diameters=core_diameters.magnitude,
        shell_width=shell_width.magnitude,
        measure_string=measure_string
    )

    # Get data from PyMieSim
    pymiesim_data = get_pymiesim_data(
        source=source,
        core_index=core_index,
        shell_index=shell_index,
        core_diameters=core_diameters,
        shell_width=shell_width,
        measure_string=measure_string
    )

    return pymiesim_data, pymiescatt_data


@pytest.mark.parametrize('measure_string', ['Qext', 'Qsca', 'Qabs', 'g', 'Qpr'])
def test_comparison(measure_string):
    """
    Test comparison between PyMieSim and PyMieScatt data for various measurements.

    Parameters:
        measure_string (str): The type of measurement to compare.
    """
    pymiesim_data, pymiescatt_data = get_comparison(
        wavelength=632 * nanometer,  # Wavelength in meters (e.g., 632 nm)
        core_index=1.4 + 0.3j * RIU,  # Complex refractive index for core
        shell_index=1.3 * RIU,  # Refractive index for shell
        core_diameters=np.geomspace(10, 6000, 800) * nanometer,  # Log-spaced array of core diameters
        shell_width=600 * nanometer,  # Shell thickness
        measure_string=measure_string
    )

    # Compare the data with a relative tolerance of 0.0001 (0.01%)
    discrepancy = np.isclose(pymiesim_data, pymiescatt_data, atol=0, rtol=1e-4)

    # Ensure that more than half of the points are within the acceptable tolerance
    assert discrepancy.astype(int).mean() > 0.5, f"Mismatch in PyMieSim vs PyMieScatt for {measure_string}"


if __name__ == "__main__":
    pytest.main([__file__])
