#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pytest
from PyMieSim import single, experiment
from PyMieSim.experiment import measure


@pytest.fixture
def source_single():
    """
    Fixture to create a Gaussian light source for single scattering simulation.

    Returns:
        single.source.Gaussian: Gaussian light source for single scattering.
    """
    return single.source.Gaussian(
        wavelength=1e-6,   # Wavelength in meters (e.g., 1 micron)
        polarization=0,    # Polarization angle
        optical_power=1,   # Optical power in watts
        NA=0.3             # Numerical aperture
    )


@pytest.fixture
def source_experiment():
    """
    Fixture to create a Gaussian light source for experiment-based simulation.

    Returns:
        experiment.source.Gaussian: Gaussian light source for experiment-based scattering.
    """
    return experiment.source.Gaussian(
        wavelength=1e-6,
        polarization=0,
        optical_power=1,
        NA=0.3
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
        diameter=1e-6,         # Diameter of the sphere in meters
        index=1.5 + 0.5j,      # Complex refractive index of the scatterer
        source=source_single,  # Associated light source
        medium_index=1         # Refractive index of the medium (e.g., air)
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
        diameter=1e-6,
        index=1.5 + 0.5j,
        source=source_experiment,
        medium_index=1
    )


def test_detector_single_polarization_filter(source_single, scatterer_single):
    """
    Test the effect of a 0-degree and 180-degree polarization filter on a photodiode detector.

    Parameters:
        source_single (single.source.Gaussian): Gaussian light source for single scattering.
        scatterer_single (single.scatterer.Sphere): Spherical scatterer.
    """
    # Create two photodiode detectors with different polarization filters (0° and 180°)
    detector_0 = single.detector.Photodiode(
        NA=0.1,
        gamma_offset=0,
        phi_offset=90,
        polarization_filter=0
    )

    detector_180 = single.detector.Photodiode(
        NA=0.1,
        gamma_offset=0,
        phi_offset=90,
        polarization_filter=180
    )

    # Calculate the coupling for both detectors
    coupling_0 = detector_0.coupling(scatterer_single)
    coupling_180 = detector_180.coupling(scatterer_single)

    # Assert that the coupling values for 0° and 180° polarization are nearly equal
    assert np.isclose(coupling_0, coupling_180, atol=1e-5), (
        f'Mismatch in coupling values for 0° and 180° polarization: {coupling_0} vs {coupling_180}'
    )


def test_detector_single_rotation(source_single, scatterer_single):
    """
    Test the effect of 0-degree and 180-degree rotation on a coherent mode detector.

    Parameters:
        source_single (single.source.Gaussian): Gaussian light source for single scattering.
        scatterer_single (single.scatterer.Sphere): Spherical scatterer.
    """
    # Create two coherent mode detectors with different rotation angles (0° and 180°)
    detector_0 = single.detector.CoherentMode(
        mode_number='LP11',
        NA=0.1,
        gamma_offset=0,
        phi_offset=40,
        polarization_filter=None,
        rotation=0
    )

    detector_180 = single.detector.CoherentMode(
        mode_number='LP11',
        NA=0.1,
        gamma_offset=0,
        phi_offset=40,
        polarization_filter=None,
        rotation=180
    )

    # Calculate the coupling for both detectors
    coupling_0 = detector_0.coupling(scatterer_single)
    coupling_180 = detector_180.coupling(scatterer_single)

    # Assert that the coupling values for 0° and 180° rotation are nearly equal
    assert np.isclose(coupling_0, coupling_180, atol=1e-5), (
        f'Mismatch in coupling values for 0° and 180° rotation: {coupling_0} vs {coupling_180}'
    )


def test_detector_experiment_polarization_filter(source_experiment, scatterer_experiment):
    """
    Test the effect of a polarization filter array on an experiment-based photodiode detector.

    Parameters:
        source_experiment (experiment.source.Gaussian): Gaussian light source for experiment-based scattering.
        scatterer_experiment (experiment.scatterer.Sphere): Spherical scatterer.
    """
    # Create a photodiode detector with two polarization filters (0° and 180°)
    detector = experiment.detector.Photodiode(
        NA=0.1,
        gamma_offset=0,
        phi_offset=90,
        polarization_filter=[0, 180],  # List of polarization filters
        sampling=100  # Sampling points for the simulation
    )

    # Setup the experiment with the scatterer, detector, and source
    setup = experiment.Setup(
        scatterer=scatterer_experiment,
        detector=detector,
        source=source_experiment
    )

    # Get the coupling values for both polarization filters
    coupling_values = setup.get(measure=measure.coupling, export_as_numpy=True).squeeze()

    # Assert that the coupling values for 0° and 180° polarization are equal
    assert coupling_values[0] == coupling_values[-1], (
        'Mismatch in coupling values for 0° and 180° polarization filters.'
    )


if __name__ == "__main__":
    pytest.main([__file__])
