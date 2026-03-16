#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pytest
from PyMieSim.units import ureg

from PyMieSim import single, experiment
from PyMieSim.polarization import PolarizationState


@pytest.fixture
def source_experiment():
    """
    Fixture to create a Gaussian light source for experiment-based simulation.

    Returns:
        experiment.source.Gaussian: Gaussian light source for experiment-based scattering.
    """
    return experiment.source_set.GaussianSet(
        wavelength=[1000] * ureg.nanometer,
        polarization=experiment.polarization_set.PolarizationSet(angles=0 * ureg.degree),
        optical_power=[1] * ureg.watt,
        numerical_aperture=[0.3],
    )


@pytest.fixture
def scatterer_single():
    """
    Fixture to create a spherical scatterer for single scattering simulation.

    Parameters:
        source_single (single.source.Gaussian): Light source for single scattering.

    Returns:
        single.scatterer.Sphere: Scatterer object for single scattering.
    """
    return single.scatterer.Sphere(
        diameter=1000 * ureg.nanometer,
        material =(1.5 + 0.5j),
        medium=1,
    )


@pytest.fixture
def scatterer_experiment():
    """
    Fixture to create a spherical scatterer for experiment-based scattering.

    Parameters:
        source_experiment (experiment.source.GaussianSet): Light source for the experiment.

    Returns:
        experiment.scatterer.SphereSet: Scatterer object for experiment-based scattering.
    """
    return experiment.scatterer_set.SphereSet(
        diameter=[1000] * ureg.nanometer,
        material=[(1.5 + 0.5j)],
        medium=[1],
    )


def test_detector_single_polarization_filter():
    """
    Test the effect of a 0-ureg.degree and 180-ureg.degree polarization filter on a photodiode detector.

    Parameters:
        source_single (single.source.Gaussian): Gaussian light source for single scattering.
        scatterer_single (single.scatterer.Sphere): Spherical scatterer.
    """
    scatterer_single = single.scatterer.Sphere(
        diameter=1000 * ureg.nanometer,
        material =(1.5 + 0.5j),
        medium=1,
    )

    source_single = single.source.Gaussian(
        wavelength=1000 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3,
    )

    # Create two photodiode detectors with different polarization filters (0° and 180°)
    detector_0 = single.detector.Photodiode(
        numerical_aperture=0.1,
        gamma_offset=0 * ureg.degree,
        phi_offset=90 * ureg.degree,
        polarization_filter=0 * ureg.degree,
        medium=1.0
    )

    detector_180 = single.detector.Photodiode(
        numerical_aperture=0.1,
        gamma_offset=0 * ureg.degree,
        phi_offset=90 * ureg.degree,
        polarization_filter=180 * ureg.degree,
        medium=1.0
    )

    setup = single.Setup(
        scatterer=scatterer_single,
        source=source_single,
        detector=detector_0
    )

    # Calculate the coupling for both detectors
    coupling_0 = setup.get("coupling")

    setup = single.Setup(
        scatterer=scatterer_single,
        source=source_single,
        detector=detector_180
    )

    coupling_180 = setup.get("coupling")

    # Assert that the coupling values for 0° and 180° polarization are nearly equal
    assert np.isclose(
        coupling_0, coupling_180, atol=1e-5
    ), f"Mismatch in coupling values for 0° and 180° polarization: {coupling_0} vs {coupling_180}"


def test_detector_single_rotation():
    """
    Test the effect of 0-ureg.degree and 180-ureg.degree rotation on a coherent mode detector.

    Parameters:
        source_single (single.source.Gaussian): Gaussian light source for single scattering.
        scatterer_single (single.scatterer.Sphere): Spherical scatterer.
    """
    # Create two coherent mode detectors with different rotation angles (0° and 180°)
    scatterer_single = single.scatterer.Sphere(
        diameter=1000 * ureg.nanometer,
        material =(1.5 + 0.5j),
        medium=1,
    )

    source_single = single.source.Gaussian(
        wavelength=1000 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3,
    )

    rotation = 0 * ureg.degree
    detector_0 = single.detector.CoherentMode(
        mode_number="LP11",
        numerical_aperture=0.1,
        gamma_offset=0 * ureg.degree,
        phi_offset=40 * ureg.degree,
        rotation=rotation,
        sampling=3000
    )
    detector_180 = single.detector.CoherentMode(
        mode_number="LP11",
        numerical_aperture=0.1,
        gamma_offset=0 * ureg.degree,
        phi_offset=40 * ureg.degree,
        rotation=rotation + 180 * ureg.degree,
        sampling=3000
    )

    setup = single.Setup(
        scatterer=scatterer_single,
        source=source_single,
        detector=detector_0
    )

    # Calculate the coupling for both detectors
    coupling_0 = setup.get("coupling")

    setup = single.Setup(
        scatterer=scatterer_single,
        source=source_single,
        detector=detector_180
    )

    coupling_180 = setup.get("coupling")

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
    detector = experiment.detector_set.PhotodiodeSet(
        numerical_aperture=[0.1],
        gamma_offset=[0] * ureg.degree,
        phi_offset=[90] * ureg.degree,
        polarization_filter=[0, 180] * ureg.degree,  # List of polarization filters
        sampling=[100],  # Sampling points for the simulation
    )

    # Setup the experiment with the scatterer, detector, and source
    setup = experiment.Setup(
        scatterer_set=scatterer_experiment, detector_set=detector, source_set=source_experiment
    )

    # Get the coupling values for both polarization filters
    dataframe = setup.get("coupling")

    # Assert that the coupling values for 0° and 180° polarization are equal
    assert (
        dataframe.coupling.values[0] == dataframe.coupling.values[-1]
    ), "Mismatch in coupling values for 0° and 180° polarization filters."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
