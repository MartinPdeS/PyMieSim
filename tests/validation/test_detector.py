#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pytest
from PyMieSim import single
from PyMieSim import experiment
from PyMieSim.experiment import measure


@pytest.fixture
def source_single():
    """Fixture for creating a Gaussian source reused across tests."""
    return single.source.Gaussian(
        wavelength=1e-6,
        polarization=0,
        optical_power=1,
        NA=0.3
    )


@pytest.fixture
def source_experiment():
    """Fixture for creating a Gaussian source reused across tests."""
    return experiment.source.Gaussian(
        wavelength=1e-6,
        polarization=0,
        optical_power=1,
        NA=0.3
    )


@pytest.fixture
def scatterer_single(source_single):
    """Fixture for creating a Gaussian source reused across tests."""
    return single.scatterer.Sphere(
        diameter=1e-6,
        index=1.5 + 0.5j,
        source=source_single,
        medium_index=1
    )


@pytest.fixture
def scatterer_experiment(source_experiment):
    """Fixture for creating a Gaussian source reused across tests."""
    return experiment.scatterer.Sphere(
        diameter=1e-6,
        index=1.5 + 0.5j,
        source=source_experiment,
        medium_index=1
    )


def test_detector_single_polarization_filter(source_single, scatterer_single):
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

    coupling_0 = detector_0.coupling(scatterer_single)
    coupling_180 = detector_180.coupling(scatterer_single)
    if not numpy.isclose(coupling_0, coupling_180, atol=1e-5):
        raise ValueError(f'Mismatch with coupling value for detector with rotation of 0 and 180 degrees. [{coupling_0} vs {coupling_180}]')


def test_detector_single_rotation(source_single, scatterer_single):
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

    coupling_0 = detector_0.coupling(scatterer_single)
    coupling_180 = detector_180.coupling(scatterer_single)
    if not numpy.isclose(coupling_0, coupling_180, atol=1e-5):
        raise ValueError(f'Mismatch with coupling value for detector with rotation of 0 and 180 degrees. [{coupling_0} vs {coupling_180}]')


def test_detector_experiment_polarization_filter(source_experiment, scatterer_experiment):
    detector = experiment.detector.Photodiode(
        NA=0.1,
        gamma_offset=0,
        phi_offset=90,
        polarization_filter=[0, 180],
        sampling=100
    )

    setup = experiment.Setup(
        scatterer=scatterer_experiment,
        detector=detector,
        source=source_experiment
    )

    coupling_values = setup.get(measure=measure.coupling, export_as_numpy=True).squeeze()

    if not coupling_values[0] == coupling_values[-1]:
        raise ValueError('Mismatch with coupling value for detector with polarization filter of 0 and 180 degrees')


if __name__ == "__main__":
    pytest.main()

# -
