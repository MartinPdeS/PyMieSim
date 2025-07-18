#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyOptik import Material
from unittest.mock import patch
import matplotlib.pyplot as plt
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyMieSim.single import plot_system

# Core, shell, and medium parameters
core_property = [Material.iron, Material.BK7, 1.0 * RIU]
shell_property = [Material.BK7, Material.iron, 1.7 * RIU]
medium_property = [Material.water, 1.4 * RIU]

attributes = [
    "size_parameter", "radius", "volume", "cross_section", "g", "Qsca", "Qext", "Qabs", "Qback", "Qratio", "an", "bn",
    "Qforward", "Qpr", "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cforward", "Cpr"
]

plotting_functions = ["get_far_field", "get_stokes", "get_spf", "get_s1s2"]


# Reusable fixture for Gaussian source
@pytest.fixture()
def source():
    return Gaussian(
        wavelength=750 * nanometer,
        polarization=0 * degree,
        optical_power=1 * watt,
        NA=0.3 * AU
    )

# Parametrized test for coupling
@pytest.mark.parametrize('core_material', core_property, ids=[f'Core:{m}' for m in core_property])
@pytest.mark.parametrize('shell_material', shell_property, ids=[f'Shell:{m}' for m in shell_property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
def test_coupling(core_material, shell_material, medium_property, source):
    detector = Photodiode(NA=0.2 * AU, gamma_offset=0 * degree, phi_offset=0 * degree)
    scatterer = CoreShell(
        core_diameter=100 * nanometer,
        shell_thickness=200 * nanometer,
        source=source,
        medium_property=medium_property,
        core_property=core_material,
        shell_property=shell_material
    )

    # Calculate optical coupling
    coupling = detector.get_coupling(scatterer)
    assert coupling is not None


# Parametrized test for attributes
@pytest.mark.parametrize('core_material', core_property, ids=[f'Core:{m}' for m in core_property])
@pytest.mark.parametrize('shell_material', shell_property, ids=[f'Shell:{m}' for m in shell_property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
@pytest.mark.parametrize('attribute', attributes)
def test_attribute(attribute, core_material, shell_material, medium_property, source):
    scatterer = CoreShell(
        core_diameter=100 * nanometer,
        shell_thickness=200 * nanometer,
        source=source,
        medium_property=medium_property,
        core_property=core_material,
        shell_property=shell_material
    )

    # Access and verify the attribute
    attr_value = getattr(scatterer, attribute)
    assert attr_value is not None


# Parametrized test for plotting functions
@pytest.mark.parametrize('core_material', core_property, ids=[f'Core:{m}' for m in core_property])
@pytest.mark.parametrize('shell_material', shell_property, ids=[f'Shell:{m}' for m in shell_property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
@pytest.mark.parametrize('plotting_function', plotting_functions)
@patch('pyvista.Plotter.show')
@patch('matplotlib.pyplot.show')
def test_plottings(mock_show_plt, mock_show_pyvista, plotting_function, core_material, shell_material, medium_property, source):
    scatterer = CoreShell(
        core_diameter=100 * nanometer,
        shell_thickness=200 * nanometer,
        source=source,
        medium_property=medium_property,
        core_property=core_material,
        shell_property=shell_material
    )

    data = getattr(scatterer, plotting_function)()
    assert data is not None
    data.plot()
    plt.close()


@patch('pyvista.Plotter.show')
def test_plot_system(mock_show, source):
    scatterer = CoreShell(
        core_diameter=100 * nanometer,
        shell_thickness=200 * nanometer,
        source=source,
        medium_property=1.0 * RIU,
        core_property=1.4 * RIU,
        shell_property=1.5 * RIU
    )

    plot_system(scatterer)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
