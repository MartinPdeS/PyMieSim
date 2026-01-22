#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode

refractive_index = [1.5 * ureg.RIU, 1.6 * ureg.RIU]
medium_refractive_index = [1.3 * ureg.RIU, 1.4 * ureg.RIU]


# Attributes to check
attributes = [
    "a1n",
    "a2n",
    "b1n",
    "b2n",
    "size_parameter",
    "cross_section",
    "g",
    "Qsca",
    "Qext",
    "Qabs",
    "Csca",
    "Cext",
    "Cabs",
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
def test_cylinder_coupling(refractive_index, medium_refractive_index, source):
    detector = Photodiode(
        NA=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        medium_refractive_index=1.0 * ureg.RIU
    )

    source = Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1 * ureg.watt,
        NA=0.3 * ureg.AU,
    )

    scatterer = InfiniteCylinder(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_refractive_index=medium_refractive_index,
        refractive_index=refractive_index,
    )

    # Calculate optical coupling
    _ = detector.get_coupling(scatterer)


@pytest.mark.parametrize("refractive_index", refractive_index, ids=[f"refractive_index:{m}" for m in refractive_index])
@pytest.mark.parametrize(
    "medium_refractive_index", medium_refractive_index, ids=[f"Medium:{m}" for m in medium_refractive_index]
)
@pytest.mark.parametrize("attribute", attributes)
def test_cylinder_attributes(attribute, refractive_index, medium_refractive_index, source):
    scatterer = InfiniteCylinder(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_refractive_index=medium_refractive_index,
        refractive_index=refractive_index,
    )

    # Access and validate the specified attribute
    _ = getattr(scatterer, attribute)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
