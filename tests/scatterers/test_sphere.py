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

methods = [
    "an",
    "bn",
]

attributes = [
    "size_parameter",
    "area",
    "Qback",
    "Qratio",
    "Qpr",
    "Csca",
    "Cext",
    "Cabs",
    "Cpr",
    "Cback",
]

plottings = [
    "get_far_field",
    "get_stokes",
    "get_spf",
    "get_s1s2",
]


@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('method', methods)
def test_sphere_method(method: str, core_type: object):
    source = PlaneWave(
        wavelength=1e-6,
        linear_polarization=0,
        amplitude=1
    )

    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        **core_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, method)()


@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('attribute', attributes)
def test_sphere_attribute(attribute: str, core_type: object):
    source = PlaneWave(
        wavelength=1e-6,
        linear_polarization=0,
        amplitude=1
    )

    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        **core_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, attribute)


@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('plotting', plottings)
def test_sphere_plottings(plotting: str, core_type: object):
    source = PlaneWave(
        wavelength=1e-6,
        linear_polarization=0,
        amplitude=1
    )

    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        **core_type,
        n_medium=1.0
    )

    data = getattr(scatterer, plotting)()

    data.plot()


# -
