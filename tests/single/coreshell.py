#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode

# Core, shell, and medium parameters
core_refractive_index = [1.2 * ureg.RIU, 1.4 * ureg.RIU, 1.0 * ureg.RIU]
shell_refractive_index = [1.8 * ureg.RIU, 1.1 * ureg.RIU, 1.7 * ureg.RIU]
medium_refractive_index = [1.3 * ureg.RIU, 1.4 * ureg.RIU]

attributes = [
    "size_parameter",
    "radius",
    "volume",
    "cross_section",
    "g",
    "Qsca",
    "Qext",
    "Qabs",
    "Qback",
    "Qratio",
    "an",
    "bn",
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


# Reusable fixture for Gaussian source
@pytest.fixture()
def source():
    return Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1 * ureg.watt,
        NA=0.3 * ureg.AU,
    )


# Parametrized test for coupling
@pytest.mark.parametrize(
    "core_material", core_refractive_index, ids=[f"Core:{m}" for m in core_refractive_index]
)
@pytest.mark.parametrize(
    "shell_material", shell_refractive_index, ids=[f"Shell:{m}" for m in shell_refractive_index]
)
@pytest.mark.parametrize(
    "medium_refractive_index", medium_refractive_index, ids=[f"Medium:{m}" for m in medium_refractive_index]
)
def test_coupling(core_material, shell_material, medium_refractive_index, source):
    detector = Photodiode(
        NA=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        medium_refractive_index=1.0 * ureg.RIU
    )
    scatterer = CoreShell(
        core_diameter=100 * ureg.nanometer,
        shell_thickness=200 * ureg.nanometer,
        source=source,
        medium_refractive_index=medium_refractive_index,
        core_refractive_index=core_material,
        shell_refractive_index=shell_material,
    )

    # Calculate optical coupling
    coupling = detector.get_coupling(scatterer)
    assert coupling is not None


# Parametrized test for attributes
@pytest.mark.parametrize(
    "core_material", core_refractive_index, ids=[f"Core:{m}" for m in core_refractive_index]
)
@pytest.mark.parametrize(
    "shell_material", shell_refractive_index, ids=[f"Shell:{m}" for m in shell_refractive_index]
)
@pytest.mark.parametrize(
    "medium_refractive_index", medium_refractive_index, ids=[f"Medium:{m}" for m in medium_refractive_index]
)
@pytest.mark.parametrize("attribute", attributes)
def test_attribute(attribute, core_material, shell_material, medium_refractive_index, source):
    scatterer = CoreShell(
        core_diameter=100 * ureg.nanometer,
        shell_thickness=200 * ureg.nanometer,
        source=source,
        medium_refractive_index=medium_refractive_index,
        core_refractive_index=core_material,
        shell_refractive_index=shell_material,
    )

    # Access and verify the attribute
    attr_value = getattr(scatterer, attribute)
    assert attr_value is not None


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
