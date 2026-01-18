#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode

# Define the core configurations for testing, now separated 'id' for clarity in tests
refractive_index = [1.2 * ureg.RIU, 1.6 * ureg.RIU]
medium_refractive_index = [1.1 * ureg.RIU, 1.4 * ureg.RIU]

attributes = [
    "size_parameter",
    "cross_section",
    "g",
    "an",
    "bn",
    "Qsca",
    "Qext",
    "Qabs",
    "Qback",
    "Qratio",
    "Qforward",
    "Qpr",
    "Csca",
    "Cext",
    "Cabs",
    "Cback",
    "Cratio",
    "Cforward",
    "Cpr",
]


@pytest.fixture()
def source():
    return Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1 * ureg.watt,
        NA=0.3 * ureg.AU,
    )


@pytest.mark.parametrize("refractive_index", refractive_index, ids=[f"refractive_index:{m}" for m in refractive_index])
@pytest.mark.parametrize(
    "medium_refractive_index", medium_refractive_index, ids=[f"Medium:{m}" for m in medium_refractive_index]
)
@pytest.mark.parametrize("attribute", attributes)
def test_sphere_attribute(attribute, refractive_index, medium_refractive_index, source):
    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_refractive_index=medium_refractive_index,
        refractive_index=refractive_index,
    )
    _ = getattr(scatterer, attribute)

    # scatterer.print_properties()


@pytest.mark.parametrize("refractive_index", refractive_index, ids=[f"refractive_index:{m}" for m in refractive_index])
@pytest.mark.parametrize(
    "medium_refractive_index", medium_refractive_index, ids=[f"Medium:{m}" for m in medium_refractive_index]
)
def test_sphere_coupling(
    refractive_index,
    medium_refractive_index,
    source,
):
    detector = Photodiode(
        sampling=100,
        NA=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        medium_refractive_index=1.0 * ureg.RIU
    )

    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_refractive_index=medium_refractive_index,
        refractive_index=refractive_index,
    )
    _ = detector.get_coupling(scatterer)


@pytest.mark.parametrize("refractive_index", refractive_index, ids=[f"refractive_index:{m}" for m in refractive_index])
@pytest.mark.parametrize(
    "medium_refractive_index", medium_refractive_index, ids=[f"Medium:{m}" for m in medium_refractive_index]
)
def test_unstructured_array_functions(refractive_index, medium_refractive_index, source):
    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_refractive_index=medium_refractive_index,
        refractive_index=refractive_index,
    )

    phi = numpy.linspace(0, numpy.pi, 5) * ureg.radian
    s1, s2 = scatterer.get_s1s2(phi)
    assert s1.shape == phi.shape
    assert s2.shape == phi.shape

    theta = numpy.linspace(0, numpy.pi / 2, 5) * ureg.radian
    I, Q, U, V = scatterer.get_stokes_parameters(phi, theta, 1 * ureg.meter)
    assert I.shape == phi.shape
    assert Q.shape == phi.shape
    assert U.shape == phi.shape
    assert V.shape == phi.shape

    E_para, E_perp = scatterer.get_farfields(phi, theta, 1 * ureg.meter)
    assert E_para.shape == phi.shape
    assert E_perp.shape == phi.shape


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
