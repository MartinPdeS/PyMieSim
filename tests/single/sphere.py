#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import Setup
from PyMieSim.material import SellmeierMaterial, SellmeierMedium

bk7 = SellmeierMaterial("BK7")

water = SellmeierMedium("water")

materials = [1.2 * ureg.RIU, bk7]
mediums = [1.1 * ureg.RIU, water]


@pytest.fixture()
def source():
    return Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )


@pytest.mark.parametrize("material", materials, ids=[f"material:{m}" for m in materials])
@pytest.mark.parametrize(
    "medium", mediums, ids=[f"Medium:{m}" for m in mediums]
)
@pytest.mark.parametrize("attribute", Sphere.property_names + ["coupling"])
def test_sphere_attribute(attribute, material, medium, source):
    detector = Photodiode(
        sampling=100,
        numerical_aperture=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        medium=1.0 * ureg.RIU
    )

    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        medium=medium,
        material=material,
    )

    setup = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector,
    )

    _ = setup.get(attribute)


@pytest.mark.parametrize("material", materials, ids=[f"material:{m}" for m in materials])
@pytest.mark.parametrize(
    "medium", mediums, ids=[f"Medium:{m}" for m in mediums]
)
def test_unstructured_array_functions(material, medium, source):
    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        medium=medium,
        material=material,
    )

    setup = Setup(
        scatterer=scatterer,
        source=source,
    )

    phi = numpy.linspace(0, numpy.pi, 5) * ureg.radian
    s1, s2 = setup.get_s1s2(phi)
    assert s1.shape == phi.shape
    assert s2.shape == phi.shape

    theta = numpy.linspace(0, numpy.pi / 2, 5) * ureg.radian
    I, Q, U, V = setup.get_stokes(phi, theta, 1 * ureg.meter)
    assert I.shape == phi.shape
    assert Q.shape == phi.shape
    assert U.shape == phi.shape
    assert V.shape == phi.shape

    E_para, E_perp = setup.get_farfields(phi, theta, 1 * ureg.meter)
    assert E_para.shape == phi.shape
    assert E_perp.shape == phi.shape


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
