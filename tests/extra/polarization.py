#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment import Setup
from PyMieSim.experiment.source import Gaussian, PolarizationSet



def test_init_linear():
    """
    Test initialization of Linear polarization.
    """
    output = PolarizationSet(angles=[50, 20] * ureg.degree)
    assert output is not None, "Initialization of Linear polarization failed!"


def test_fail_init_linear():
    """
    Test fail of initialization of Linear polarization.
    """
    with pytest.raises(AttributeError):
        PolarizationSet(angles=[50, 20])


def test_init_jones_vector():
    """
    Test initialization of Jones Vector polarization.
    """
    output = PolarizationSet(jones_vectors=[(1, 0), (0, 1)])
    assert output is not None, "Initialization of Jones Vector polarization failed!"


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
