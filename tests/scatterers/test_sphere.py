#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyOptik import Material
from unittest.mock import patch
import matplotlib.pyplot as plt

# Define the core configurations for testing, now separated 'id' for clarity in tests
core_configs = [
    {'config': {'material': Material.BK7}, 'id': 'core:BK7'},
    {'config': {'index': 1.6}, 'id': 'core:1.6'}
]

medium_configs = [
    {'config': {'medium_material': Material.water}, 'id': 'medium:water'},
    {'config': {'medium_index': 1.4}, 'id': 'medium:index'}
]

methods = ['an', 'bn']

attributes = [
    "size_parameter", "area", "g",
    "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qforward", "Qpr",
    "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cforward", "Cpr"
]

plotting_functions = [
    "get_far_field", "get_stokes", "get_spf", "get_s1s2",
]


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('method', methods)
def test_sphere_method(method, core_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )
    # Pass only the actual configuration dictionary to the constructor
    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        **medium_config['config'],
        **core_config['config']
    )

    _ = getattr(scatterer, method)()

    _ = getattr(scatterer, method)(max_order=3)


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('attribute', attributes)
def test_sphere_attribute(attribute, core_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )
    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        **medium_config['config'],
        **core_config['config']
    )
    _ = getattr(scatterer, attribute)

    scatterer.print_properties()


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
def test_sphere_coupling(core_config, medium_config):
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
    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        **medium_config['config'],
        **core_config['config']
    )
    _ = detector.coupling(scatterer)


@pytest.mark.parametrize('core_config', core_configs, ids=[config['id'] for config in core_configs])
@pytest.mark.parametrize('medium_config', medium_configs, ids=[config['id'] for config in medium_configs])
@pytest.mark.parametrize('plotting_function', plotting_functions)
@patch('pyvista.Plotter.show')
@patch('matplotlib.pyplot.show')
def test_sphere_plottings(mock_show_plt, mock_show_pyvista, plotting_function, core_config, medium_config):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )
    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        **medium_config['config'],
        **core_config['config']
    )
    data = getattr(scatterer, plotting_function)()
    data.plot()
    assert data is not None, "Plotting data should not be None"

    plt.close()


if __name__ == "__main__":
    pytest.main([__file__])


# -
