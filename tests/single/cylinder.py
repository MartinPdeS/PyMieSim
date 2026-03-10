#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import PolarizationState
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import Setup

materials = [1.5 * ureg.RIU, 1.6 * ureg.RIU]
mediums = [1.3 * ureg.RIU, 1.4 * ureg.RIU]

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
@pytest.mark.parametrize("attribute", InfiniteCylinder.property_names + ["coupling"])
def test_cylinder_coupling(attribute, material, medium, source):
    detector = Photodiode(
        numerical_aperture=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        medium=medium
    )

    source = Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )

    scatterer = InfiniteCylinder(
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



if __name__ == "__main__":
    pytest.main(["-W error", __file__])
