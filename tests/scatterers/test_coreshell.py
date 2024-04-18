#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.materials import Silver, BK7, Aluminium

# Core and shell configurations with clear separation of test ids and parameters
core_configs = [
    {'config': {'core_material': BK7}, 'id': 'BK7'},
    {'config': {'core_material': Silver}, 'id': 'Silver'},
    {'config': {'core_index': 1.4}, 'id': 'Index 1.4'}
]

shell_configs = [
    {'config': {'shell_material': BK7}, 'id': 'BK7'},
    {'config': {'shell_material': Silver}, 'id': 'Silver'},
    {'config': {'shell_material': Aluminium}, 'id': 'Aluminium'},
    {'config': {'shell_index': 1.4}, 'id': 'Index 1.4'}
]

# Methods and attributes to test
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


@pytest.mark.parametrize('shell_config', shell_configs, ids=[config['id'] for config in shell_configs])
@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('method', methods)
def test_coreshell_method(method, core_config, shell_config):
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
        medium_index=1.0,
        **core_config['config'],
        **shell_config['config']
    )

    # Execute the method
    _ = getattr(scatterer, method)()


@pytest.mark.parametrize('shell_config', shell_configs, ids=[config['id'] for config in shell_configs])
@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
def test_coreshell_coupling(core_config, shell_config):
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
    scatterer = CoreShell(
        core_diameter=100e-9,
        shell_width=200e-9,
        source=source,
        medium_index=1.0,
        **core_config['config'],
        **shell_config['config']
    )

    # Calculate optical coupling
    _ = detector.coupling(scatterer)


@pytest.mark.parametrize('shell_config', shell_configs, ids=[config['id'] for config in shell_configs])
@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('attribute', attributes)
def test_coreshell_attribute(attribute, core_config, shell_config):
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
        medium_index=1.0,
        **core_config['config'],
        **shell_config['config']
    )

    # Access and verify the attribute
    _ = getattr(scatterer, attribute)

# -
