#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.scatterer import Cylinder
from PyMieSim.source import PlaneWave
from PyMieSim.materials import Silver, BK7

cores_type = [
    {'material': BK7},
    {'material': Silver},
    {'index': 1.4}
]

methods = [
    "get_far_field",
    "get_stokes",
    "get_spf",
    "get_s1s2",
    "a1n",
    "a2n",
    "b1n",
    "b2n",
]


attributes = [
    "size_parameter",
    "area",
    # "Qsca",
    # "Qext",
    # "Qratio",
    # "Qpr",
    # "Csca",
    # "Cext",
    # "Cratio",
    # "Cpr",
]


@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('method', methods)
def test_cylinder_methods(method, core_type):
    source = PlaneWave(
        wavelength=1e-6,
        linear_polarization=0,
        amplitude=1
    )

    scatterer = Cylinder(
        diameter=100e-9,
        source=source,
        **core_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, method)()


@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('attribute', attributes)
def test_cylinder_attributes(attribute, core_type):
    source = PlaneWave(
        wavelength=1e-6,
        linear_polarization=0,
        amplitude=1
    )

    scatterer = Cylinder(
        diameter=100e-9,
        source=source,
        **core_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, attribute)


# -
