#!/usr/bin/env python
"""Regression tests for the high-level experiment API."""

import numpy as np
import pytest

from PyMieSim.experiment import Setup
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.units import ureg


def _setup():
    source = GaussianSet(
        wavelength=[600, 700] * ureg.nanometer,
        polarization=PolarizationSet(angles=0 * ureg.degree),
        optical_power=[1e-3] * ureg.watt,
        numerical_aperture=[0.2],
    )
    scatterer = SphereSet(
        diameter=[100, 200] * ureg.nanometer,
        material=[1.4],
        medium=[1.0],
    )
    return Setup(scatterer_set=scatterer, source_set=source)


def test_get_rejects_empty_measure_list():
    with pytest.raises(ValueError, match="At least one measure"):
        _setup().get()


def test_get_rejects_unknown_measure():
    with pytest.raises(ValueError, match="Unknown measure"):
        _setup().get("not_a_measure")


def test_get_as_numpy_preserves_measure_order_and_shape():
    experiment = _setup()

    values = experiment.get("Qext", "Qsca", as_numpy=True)

    assert values.shape == (2, 2, 2)
    np.testing.assert_allclose(values[0], experiment.get("Qext", as_numpy=True))
    np.testing.assert_allclose(values[1], experiment.get("Qsca", as_numpy=True))


def test_get_single_numpy_measure_keeps_simulation_shape():
    values = _setup().get("Qsca", as_numpy=True)

    assert values.shape == (2, 2)
