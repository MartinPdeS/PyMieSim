#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from TypedUnit import ureg
from PyOptik import Material
import matplotlib.pyplot as plt
from unittest.mock import patch

from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import plot_system

# Core configurations separated for clarity and functionality

property = [Material.BK7, 1.6 * ureg.RIU]
medium_property = [Material.water, 1.4 * ureg.RIU]


# Attributes to check
attributes = [
    "a1n", "a2n", "b1n", "b2n",
    "size_parameter", "cross_section", "g",
    "Qsca", "Qext", "Qabs",
    "Csca", "Cext", "Cabs",
]


@pytest.fixture()
def source():
    return Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1 * ureg.watt,
        NA=0.3 * ureg.AU
    )


@pytest.mark.parametrize('property', property, ids=[f'property:{m}' for m in property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
def test_cylinder_coupling(property, medium_property, source):
    detector = Photodiode(
        NA=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
    )

    scatterer = Cylinder(
        diameter=100 * ureg.nanometer,
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
        diameter=100 * ureg.nanometer,
        source=source,
        medium_property=medium_property,
        property=property
    )

    # Access and validate the specified attribute
    _ = getattr(scatterer, attribute)

    scatterer.print_properties()


@patch('pyvista.Plotter.show')
def test_plot_system(mock_show, source):
    scatterer = Cylinder(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_property=1.0 * ureg.RIU,
        property=1.5 * ureg.RIU
    )

    plot_system(scatterer)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
