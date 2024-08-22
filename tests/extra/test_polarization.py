#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment import Setup, measure
from PyMieSim.experiment.source import Gaussian
from PyMieSim.polarization import Linear, JonesVector, RightCircular, LeftCircular


def test_init_linear():
    """
    Test initialization of Linear polarization.
    """
    output = Linear(angles=[50, 20])
    assert output is not None, 'Initialization of Linear polarization failed!'


def test_init_right_circular():
    """
    Test initialization of Right Circular polarization.
    """
    output = RightCircular()
    assert output is not None, 'Initialization of Right Circular polarization failed!'


def test_init_left_circular():
    """
    Test initialization of Left Circular polarization.
    """
    output = LeftCircular()
    assert output is not None, 'Initialization of Left Circular polarization failed!'


def test_init_jones_vector():
    """
    Test initialization of Jones Vector polarization.
    """
    output = JonesVector(elements=[(1, 0), (0, 1)])
    assert output is not None, 'Initialization of Jones Vector polarization failed!'


polarizations = [
    Linear(angles=[50, 20]),
    RightCircular(),
    LeftCircular(),
    JonesVector(elements=[(1, 0), (0, 1)])
]


@pytest.mark.parametrize('polarization_0', polarizations, ids=lambda p: p.__class__.__name__)
@pytest.mark.parametrize('polarization_1', polarizations, ids=lambda p: p.__class__.__name__)
def test_addition_operator(polarization_0, polarization_1):
    """
    Test the addition operator for different polarizations.
    """
    output = polarization_0 + polarization_1
    assert output is not None, 'Addition of polarizations failed!'


@pytest.mark.parametrize('polarization_0', polarizations, ids=lambda p: p.__class__.__name__)
@pytest.mark.parametrize('polarization_1', polarizations, ids=lambda p: p.__class__.__name__)
def test_api(polarization_0, polarization_1):
    """
    Test the API integration with different polarizations.
    """
    # Setup Gaussian source
    source = Gaussian(
        wavelength=1e-6,
        polarization=polarization_0 + polarization_1,
        optical_power=1e-3,
        NA=0.2
    )

    # Configure the spherical scatterer
    scatterer = Sphere(
        diameter=1e-6,
        source=source,
        index=1.4,
        medium_index=1.0
    )

    # Set up and run the experiment
    experiment = Setup(
        scatterer=scatterer,
        source=source,
    )

    result = experiment.get(measure.Qsca)
    assert result is not None, 'Experiment setup or measurement failed!'

if __name__ == "__main__":
    pytest.main()

# -
