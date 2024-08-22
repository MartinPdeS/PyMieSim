#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyOptik import UsualMaterial

# Core and shell configurations with clear separation of test ids and parameters
core_configs = [
    {'config': {'core_material': UsualMaterial.Aluminium}, 'id': 'Shell:Aluminium'},
    {'config': {'core_material': UsualMaterial.Silver}, 'id': 'Shell:Silver'},
    {'config': {'core_index': 1.6}, 'id': 'Shell:1.6'}
]

shell_configs = [
    {'config': {'shell_material': UsualMaterial.BK7}, 'id': 'BK7'},
    {'config': {'shell_material': UsualMaterial.Aluminium}, 'id': 'Core:Aluminium'},
    {'config': {'shell_index': 1.7}, 'id': 'Shell:1.7'}
]

medium_configs = [
    {'config': {'medium_material': UsualMaterial.BK7}, 'id': 'Medium:BK7'},
    {'config': {'medium_index': 1.4}, 'id': 'Medium:1.4'}
]

# Methods and attributes to test
methods = ["an", "bn"]

attributes = [
    "size_parameter", "area", "g",
    "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qforward", "Qpr",
    "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cforward", "Cpr",
]

plottings = [
    "get_far_field", "get_stokes", "get_spf", "get_s1s2",
]


@pytest.mark.parametrize('shell_config', shell_configs, ids=[config['id'] for config in shell_configs])
@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('method', methods)
def test_coreshell_method(method, core_config, shell_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )
    scatterer = CoreShell(
        core_diameter=100e-9,
        shell_width=200e-9,
        source=source,
        **medium_config['config'],
        **core_config['config'],
        **shell_config['config']
    )

    # Execute the method
    _ = getattr(scatterer, method)()

    _ = getattr(scatterer, method)(max_order=3)


@pytest.mark.parametrize('shell_config', shell_configs, ids=[config['id'] for config in shell_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
def test_coreshell_coupling(core_config, shell_config, medium_config):
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

    scatterer = CoreShell(
        core_diameter=100e-9,
        shell_width=200e-9,
        source=source,
        **medium_config['config'],
        **core_config['config'],
        **shell_config['config']
    )

    # Calculate optical coupling
    _ = detector.coupling(scatterer)


@pytest.mark.parametrize('shell_config', shell_configs, ids=[config['id'] for config in shell_configs])
@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('attribute', attributes)
def test_coreshell_attribute(attribute, core_config, shell_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )

    scatterer = CoreShell(
        core_diameter=100e-9,
        shell_width=200e-9,
        source=source,
        **medium_config['config'],
        **core_config['config'],
        **shell_config['config']
    )

    # Access and verify the attribute
    _ = getattr(scatterer, attribute)


@pytest.mark.parametrize('shell_config', shell_configs, ids=[config['id'] for config in shell_configs])
@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('plotting', plottings)
def test_coreshell_plottings(plotting, core_config, shell_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )

    scatterer = CoreShell(
        core_diameter=100e-9,
        shell_width=200e-9,
        source=source,
        **medium_config['config'],
        **core_config['config'],
        **shell_config['config']
    )

    data = getattr(scatterer, plotting)()
    assert data is not None, "Plotting data should not be None"


if __name__ == "__main__":
    pytest.main()

# -
