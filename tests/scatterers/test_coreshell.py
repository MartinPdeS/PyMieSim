#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.materials import Silver, BK7, Aluminium

cores_type = [
    {'core_material': BK7},
    {'core_material': Silver},
    {'core_index': 1.4}
]

shells_type = [
    {'shell_material': BK7},
    {'shell_material': Silver},
    {'shell_material': Aluminium},
    {'shell_index': 1.4}
]

methods = [
    "get_far_field",
    "get_stokes",
    "get_spf",
    "get_s1s2",
    "an",
    "bn",
]

attributes = [
    "size_parameter",
    "area",
    "Qsca",
    "Qext",
    "Qback",
    "Qratio",
    "Qpr",
    "Csca",
    "Cext",
    "Cback",
    "Cratio",
    "Cpr",
]


@pytest.mark.parametrize('shell_type', shells_type, ids=['BK7', 'Silver', 'Aluminum', 'Index'])
@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('method', methods)
def test_coreshell_method(method, core_type, shell_type):
    source = Gaussian(
        wavelength=750e-9,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )
    scatterer = CoreShell(
        core_diameter=100e-9,
        shell_width=200e-9,
        source=source,
        **core_type,
        **shell_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, method)()


@pytest.mark.parametrize('shell_type', shells_type, ids=['BK7', 'Silver', 'Aluminum', 'Index'])
@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('attribute', attributes)
def test_coreshell_attribute(attribute, core_type, shell_type):
    source = Gaussian(
        wavelength=750e-9,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )

    scatterer = CoreShell(
        core_diameter=100e-9,
        shell_width=200e-9,
        source=source,
        **core_type,
        **shell_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, attribute)


# -
