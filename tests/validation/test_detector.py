#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pytest
from PyMieSim.units import ureg

from PyMieSim import single, experiment
import PyMieSim
PyMieSim.debug_mode = False

@pytest.fixture
def source_single():
    """
    Fixture to create a Gaussian light source for single scattering simulation.

    Returns:
        single.source.Gaussian: Gaussian light source for single scattering.
    """
    return single.source.Gaussian(
        wavelength=1000 * ureg.nanometer,  # Wavelength in meters (e.g., 1 micron)
        polarization=0 * ureg.degree,  # Polarization angle
        optical_power=1 * ureg.watt,  # Optical power in ureg.watts
        NA=0.3 * ureg.AU,  # Numerical aperture
    )


@pytest.fixture
def source_experiment():
    """
    Fixture to create a Gaussian light source for experiment-based simulation.

    Returns:
        experiment.source.Gaussian: Gaussian light source for experiment-based scattering.
    """
    return experiment.source.Gaussian(
        wavelength=1000 * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1 * ureg.watt,
        NA=0.3 * ureg.AU,
    )


@pytest.fixture
def scatterer_single(source_single):
    """
    Fixture to create a spherical scatterer for single scattering simulation.

    Parameters:
        source_single (single.source.Gaussian): Light source for single scattering.

    Returns:
        single.scatterer.Sphere: Scatterer object for single scattering.
    """
    return single.scatterer.Sphere(
        diameter=1000 * ureg.nanometer,  # Diameter of the sphere in meters
        refractive_index=(1.5 + 0.5j) * ureg.RIU,  # Complex refractive index of the scatterer
        source=source_single,  # Associated light source
        medium_refractive_index=1 * ureg.RIU,  # Refractive index of the medium (e.g., air)
    )


@pytest.fixture
def scatterer_experiment(source_experiment):
    """
    Fixture to create a spherical scatterer for experiment-based scattering.

    Parameters:
        source_experiment (experiment.source.Gaussian): Light source for the experiment.

    Returns:
        experiment.scatterer.Sphere: Scatterer object for experiment-based scattering.
    """
    return experiment.scatterer.Sphere(
        diameter=1000 * ureg.nanometer,
        property=(1.5 + 0.5j) * ureg.RIU,
        source=source_experiment,
        medium_property=1 * ureg.RIU,
    )


def test_detector_single_polarization_filter(scatterer_single):
    """
    Test the effect of a 0-ureg.degree and 180-ureg.degree polarization filter on a photodiode detector.

    Parameters:
        source_single (single.source.Gaussian): Gaussian light source for single scattering.
        scatterer_single (single.scatterer.Sphere): Spherical scatterer.
    """
    # Create two photodiode detectors with different polarization filters (0° and 180°)
    detector_0 = single.detector.Photodiode(
        NA=0.1 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=90 * ureg.degree,
        polarization_filter=0 * ureg.degree,
        medium_refractive_index=1.0 * ureg.RIU
    )

    detector_180 = single.detector.Photodiode(
        NA=0.1 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=90 * ureg.degree,
        polarization_filter=180 * ureg.degree,
        medium_refractive_index=1.0 * ureg.RIU
    )

    # Calculate the coupling for both detectors
    coupling_0 = detector_0.get_coupling(scatterer_single)
    coupling_180 = detector_180.get_coupling(scatterer_single)

    # Assert that the coupling values for 0° and 180° polarization are nearly equal
    assert np.isclose(
        coupling_0, coupling_180, atol=1e-5
    ), f"Mismatch in coupling values for 0° and 180° polarization: {coupling_0} vs {coupling_180}"


def test_detector_single_rotation(scatterer_single):
    """
    Test the effect of 0-ureg.degree and 180-ureg.degree rotation on a coherent mode detector.

    Parameters:
        source_single (single.source.Gaussian): Gaussian light source for single scattering.
        scatterer_single (single.scatterer.Sphere): Spherical scatterer.
    """
    # Create two coherent mode detectors with different rotation angles (0° and 180°)

    rotation = 0 * ureg.degree
    detector_0 = single.detector.CoherentMode(
        mode_number="LP11",
        NA=0.1 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=40 * ureg.degree,
        rotation=rotation,
        sampling=3000
    )
    detector_180 = single.detector.CoherentMode(
        mode_number="LP11",
        NA=0.1 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=40 * ureg.degree,
        rotation=rotation + 180 * ureg.degree,
        sampling=3000
    )

    # Calculate the coupling for both detectors
    coupling_0 = detector_0.get_coupling(scatterer_single)
    coupling_180 = detector_180.get_coupling(scatterer_single)

    # Assert that the coupling values for 0° and 180° rotation are nearly equal
    assert np.isclose(
        coupling_0, coupling_180, rtol=1e-3
    ), f"Mismatch in coupling values for 0° and 180° rotation: {coupling_0} vs {coupling_180}"


def test_detector_experiment_polarization_filter(
    source_experiment, scatterer_experiment
):
    """
    Test the effect of a polarization filter array on an experiment-based photodiode detector.

    Parameters:
        source_experiment (experiment.source.Gaussian): Gaussian light source for experiment-based scattering.
        scatterer_experiment (experiment.scatterer.Sphere): Spherical scatterer.
    """
    # Create a photodiode detector with two polarization filters (0° and 180°)
    detector = experiment.detector.Photodiode(
        NA=0.1 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=90 * ureg.degree,
        polarization_filter=[0, 180] * ureg.degree,  # List of polarization filters
        sampling=100 * ureg.AU,  # Sampling points for the simulation
    )

    # Setup the experiment with the scatterer, detector, and source
    setup = experiment.Setup(
        scatterer=scatterer_experiment, detector=detector, source=source_experiment
    )

    # Get the coupling values for both polarization filters
    dataframe = setup.get("coupling")

    # Assert that the coupling values for 0° and 180° polarization are equal
    assert (
        dataframe.coupling.values[0] == dataframe.coupling.values[-1]
    ), "Mismatch in coupling values for 0° and 180° polarization filters."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
