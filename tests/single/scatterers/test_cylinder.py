#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyOptik import Material
from unittest.mock import patch
import matplotlib.pyplot as plt
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyMieSim.single import plot_system

# Core configurations separated for clarity and functionality

property = [Material.BK7, 1.6 * RIU]
medium_property = [Material.water, 1.4 * RIU]


# Attributes to check
attributes = [
    "a1n", "a2n", "b1n", "b2n",
    "size_parameter", "cross_section", "g",
    "Qsca", "Qext", "Qabs",
    "Csca", "Cext", "Cabs",
]

plotting_functions = [
    "get_farfield", "get_stokes", "get_spf", "get_s1s2",
]


@pytest.fixture()
def source():
    return Gaussian(
        wavelength=750 * nanometer,
        polarization=0 * degree,
        optical_power=1 * watt,
        NA=0.3 * AU
    )


@pytest.mark.parametrize('property', property, ids=[f'property:{m}' for m in property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
def test_cylinder_coupling(property, medium_property, source):
    detector = Photodiode(
        NA=0.2 * AU,
        gamma_offset=0 * degree,
        phi_offset=0 * degree,
    )

    scatterer = Cylinder(
        diameter=100 * nanometer,
        source=source,
        medium_property=medium_property,
        property=property
    )

    # Calculate optical coupling
    _ = detector.get_coupling(scatterer)


@pytest.mark.parametrize('property', property, ids=[f'property:{m}' for m in property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
@pytest.mark.parametrize('attribute', attributes)
def test_cylinder_attributes(attribute, property, medium_property, source):
    scatterer = Cylinder(
        diameter=100 * nanometer,
        source=source,
        medium_property=medium_property,
        property=property
    )

    # Access and validate the specified attribute
    _ = getattr(scatterer, attribute)

    scatterer.print_properties()


@pytest.mark.parametrize('property', property, ids=[f'property:{m}' for m in property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
@pytest.mark.parametrize('plotting_function', plotting_functions)
@patch('pyvista.Plotter.show')
@patch('matplotlib.pyplot.show')
def test_cylinder_plottings(mock_show_plt, mock_show_pyvista, plotting_function, property, medium_property, source):
    scatterer = Cylinder(
        diameter=100 * nanometer,
        source=source,
        medium_property=medium_property,
        property=property
    )

    data = getattr(scatterer, plotting_function)()
    data.plot()
    assert data is not None, "Plotting data should not be None"

    plt.close()
    import pyvista

    pyvista.close_all()


@patch('pyvista.Plotter.show')
def test_plot_system(mock_show, source):
    scatterer = Cylinder(
        diameter=100 * nanometer,
        source=source,
        medium_property=1.0 * RIU,
        property=1.5 * RIU
    )

    plot_system(scatterer)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
