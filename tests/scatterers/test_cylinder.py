#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyOptik import UsualMaterial

# Core configurations separated for clarity and functionality
core_configs = [
    {'config': {'material': UsualMaterial.Silver}, 'id': 'Core:Silver'},
    {'config': {'index': 1.6}, 'id': 'Core:1.6'}
]

medium_configs = [
    {'config': {'medium_material': UsualMaterial.BK7}, 'id': 'Medium:BK7'},
    {'config': {'medium_index': 1.4}, 'id': 'Medium:1.4'}
]

# Methods to be tested
methods = [
    "a1n",
    "a2n",
    "b1n",
    "b2n",
]

# Attributes to check
attributes = [
    "size_parameter", "area", "g",
    "Qsca", "Qext", "Qabs", "Qpr",
    "Csca", "Cext", "Cabs", "Cpr"
]

plottings = [
    "get_far_field", "get_stokes", "get_spf", "get_s1s2",
]


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('method', methods)
def test_cylinder_methods(method, core_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )
    scatterer = Cylinder(
        diameter=100e-9,
        source=source,
        **medium_config['config'],
        **core_config['config']
    )

    # Execute method
    _ = getattr(scatterer, method)()

    _ = getattr(scatterer, method)(max_order=3)


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
def test_cylinder_coupling(core_config, medium_config):
    detector = Photodiode(
        NA=0.2,
        gamma_offset=0,
        phi_offset=0,
    )

    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )

    scatterer = Cylinder(
        diameter=100e-9,
        source=source,
        **medium_config['config'],
        **core_config['config']
    )

    # Calculate optical coupling
    _ = detector.coupling(scatterer)


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('attribute', attributes)
def test_cylinder_attributes(attribute, core_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )

    scatterer = Cylinder(
        diameter=100e-9,
        source=source,
        **medium_config['config'],
        **core_config['config']
    )

    # Access and validate the specified attribute
    _ = getattr(scatterer, attribute)


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('plotting', plottings)
def test_cylinder_plottings(plotting, core_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )

    scatterer = Cylinder(
        diameter=100e-9,
        source=source,
        **medium_config['config'],
        **core_config['config']
    )

    data = getattr(scatterer, plotting)()
    assert data is not None, "Plotting data should not be None"


if __name__ == "__main__":
    pytest.main()

# -
