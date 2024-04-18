#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.materials import Silver, BK7

# Core configurations separated for clarity and functionality
core_configs = [
    {'config': {'material': BK7}, 'id': 'BK7'},
    {'config': {'material': Silver}, 'id': 'Silver'},
    {'config': {'index': 1.4}, 'id': 'Index 1.4'}
]

# Methods to be tested
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

# Attributes to check
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


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('method', methods)
def test_cylinder_methods(method, core_config):
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
        medium_index=1.0,
        **core_config['config']
    )

    # Execute method
    _ = getattr(scatterer, method)()


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
def test_cylinder_coupling(core_config):
    detector = Photodiode(
        NA=0.2,
        gamma_offset=0,
        phi_offset=0,
    )
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
        medium_index=1.0,
        **core_config['config']
    )

    # Calculate optical coupling
    _ = detector.coupling(scatterer)


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('attribute', attributes)
def test_cylinder_attributes(attribute, core_config):
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
        medium_index=1.0,
        **core_config['config']
    )

    # Access and validate the specified attribute
    _ = getattr(scatterer, attribute)

# -
