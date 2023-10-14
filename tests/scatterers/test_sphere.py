#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.scatterer import Sphere
from PyMieSim.source import PlaneWave
from PyMieSim.materials import Silver, BK7

cores_type = [
    {'material': BK7},
    {'material': Silver},
    {'index': 1.4}
]

attributes = [
    "get_far_field",
    "get_stokes",
    "get_spf",
    "get_s1s2"
]


@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('attribute', attributes)
def test_coreshell_experiment(attribute, core_type):
    source = PlaneWave(
        wavelength=1e-6,
        polarization=0,
        amplitude=1
    )

    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        **core_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, attribute)()


# -
