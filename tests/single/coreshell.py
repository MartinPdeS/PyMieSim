#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import Setup
from PyMieSim.material import SellmeierMaterial, SellmeierMedium

bk7 = SellmeierMaterial("BK7")

water = SellmeierMedium("water")

# Core, shell, and medium parameters
core_refractive_index = [1.2 * ureg.RIU, bk7, 1.0 * ureg.RIU]
shell_refractive_index = [1.8 * ureg.RIU, 1.1 * ureg.RIU, 1.7 * ureg.RIU]
medium_refractive_index = [1.3 * ureg.RIU, water]


# Reusable fixture for Gaussian source
@pytest.fixture()
def source():
    return Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )


# Parametrized test for coupling
@pytest.mark.parametrize(
    "core_material", core_refractive_index, ids=[f"Core:{m}" for m in core_refractive_index]
)
@pytest.mark.parametrize(
    "shell_material", shell_refractive_index, ids=[f"Shell:{m}" for m in shell_refractive_index]
)
@pytest.mark.parametrize(
    "medium", medium_refractive_index, ids=[f"Medium:{m}" for m in medium_refractive_index]
)
@pytest.mark.parametrize("attribute", CoreShell.property_names + ["coupling"])
def test_coupling(attribute, core_material, shell_material, medium, source):
    detector = Photodiode(
        numerical_aperture=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        medium=medium
    )
    scatterer = CoreShell(
        core_diameter=100 * ureg.nanometer,
        shell_thickness=200 * ureg.nanometer,
        medium=medium,
        core_material=core_material,
        shell_material=shell_material,
    )

    setup = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector,
    )

    # Calculate optical coupling
    coupling = setup.get(attribute)
    assert coupling is not None



if __name__ == "__main__":
    pytest.main(["-W error", __file__])
