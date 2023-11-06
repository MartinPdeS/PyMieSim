#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.scatterer import CoreShell
from PyMieSim.source import PlaneWave
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

attributes = [
    "get_far_field",
    "get_stokes",
    "get_spf",
    "get_s1s2"
]


@pytest.mark.parametrize('shell_type', shells_type, ids=['BK7', 'Silver', 'Aluminum', 'Index'])
@pytest.mark.parametrize('core_type', cores_type, ids=['BK7', 'Silver', 'Index'])
@pytest.mark.parametrize('attribute', attributes)
def test_coreshell_experiment(attribute, core_type, shell_type):
    source = PlaneWave(
        wavelength=1e-6,
        linear_polarization=0,
        amplitude=1
    )

    scatterer = CoreShell(
        core_diameter=100e-9,
        shell_width=200e-9,
        source=source,
        **core_type,
        **shell_type,
        n_medium=1.0
    )

    _ = getattr(scatterer, attribute)()


# -
