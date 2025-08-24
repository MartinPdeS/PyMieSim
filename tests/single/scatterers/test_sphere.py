#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy
from TypedUnit import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyOptik import Material
from unittest.mock import patch
import matplotlib.pyplot as plt

# Define the core configurations for testing, now separated 'id' for clarity in tests
property = [Material.BK7, 1.6 * ureg.RIU]
medium_property = [Material.water, 1.4 * ureg.RIU]

attributes = [
    "size_parameter", "cross_section", "g", 'an', 'bn',
    "Qsca", "Qext", "Qabs", "Qback", "Qratio", "Qforward", "Qpr",
    "Csca", "Cext", "Cabs", "Cback", "Cratio", "Cforward", "Cpr"
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
@pytest.mark.parametrize('attribute', attributes)
def test_sphere_attribute(attribute, property, medium_property, source):
    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_property=medium_property,
        property=property
    )
    _ = getattr(scatterer, attribute)

    scatterer.print_properties()


@pytest.mark.parametrize('property', property, ids=[f'property:{m}' for m in property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
def test_sphere_coupling(property, medium_property, source, ):
    detector = Photodiode(
        NA=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
    )

    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_property=medium_property,
        property=property
    )
    _ = detector.get_coupling(scatterer)


@pytest.mark.parametrize('property', property, ids=[f'property:{m}' for m in property])
@pytest.mark.parametrize('medium_property', medium_property, ids=[f'Medium:{m}' for m in medium_property])
def test_unstructured_array_functions(property, medium_property, source):
    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_property=medium_property,
        property=property
    )

    phi = numpy.linspace(0, numpy.pi, 5) * ureg.radian
    s1, s2 = scatterer.get_s1s2_array(phi)
    assert s1.shape == phi.shape
    assert s2.shape == phi.shape

    theta = numpy.linspace(0, numpy.pi / 2, 5) * ureg.radian
    I, Q, U, V = scatterer.get_stokes_array(phi, theta)
    assert I.shape == phi.shape
    assert Q.shape == phi.shape
    assert U.shape == phi.shape
    assert V.shape == phi.shape


    E_para, E_perp = scatterer.get_farfield_array(phi, theta)
    assert E_para.shape == phi.shape
    assert E_perp.shape == phi.shape

if __name__ == "__main__":
    pytest.main(["-W error", __file__])
