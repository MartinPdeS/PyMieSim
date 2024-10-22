#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment import Setup
from PyMieSim.experiment.source import Gaussian
from PyMieSim.polarization import Linear, JonesVector, RightCircular, LeftCircular
from PyMieSim.units import nanometer, degree, watt, AU, RIU


def test_init_linear():
    """
    Test initialization of Linear polarization.
    """
    output = Linear(element=[50, 20] * degree)
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
    output = JonesVector(element=[(1, 0), (0, 1)])
    assert output is not None, 'Initialization of Jones Vector polarization failed!'


polarizations = [
    Linear(element=[50, 20] * degree),
    RightCircular(),
    LeftCircular(),
    JonesVector(element=[(1, 0), (0, 1)])
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
    # polarization = Linear(element=[50] * degree)
    source = Gaussian(
        wavelength=np.linspace(600, 1000, 50) * nanometer,
        polarization=0 * degree,
        optical_power=1e-3 * watt,
        NA=0.2 * AU
    )

    # Configure the spherical scatterer
    scatterer = Sphere(
        diameter=np.linspace(400, 1400, 10) * nanometer,
        source=source,
        property=1.4 * RIU,
        medium_property=1.0 * RIU
    )

    # Set up and run the experiment
    experiment = Setup(scatterer=scatterer, source=source)

    result = experiment.get('Qsca')
    assert result is not None, 'Experiment setup or measurement failed!'


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
