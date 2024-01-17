#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
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
    "Qsca",
    "Qext",
    "Qabs",
    "Csca",
    "Cext",
    "Cabs",
    "g",
]


@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('method', methods)
def test_cylinder_methods(method, core_type):
    source = Gaussian(
        wavelength=750e-9,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
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
    source = Gaussian(
        wavelength=750e-9,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )

    scatterer = Cylinder(
        diameter=100e-9,
        source=source,
        **core_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, attribute)


# -
