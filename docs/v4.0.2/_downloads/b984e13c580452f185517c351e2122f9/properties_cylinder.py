"""
Print properties
================

This example demonstrates the computation of scattering properties using PyMieSim.
"""

from PyMieSim.units import ureg
from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState

polarization_state = PolarizationState(angle=30 * ureg.degree)

source = Gaussian(
    wavelength=750 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = InfiniteCylinder(
    diameter=300 * ureg.nanometer,
    material=(1.4 + 0.1j) * ureg.RIU,
    medium=1.33 * ureg.RIU,
)

scatterer.init(source)

scatterer.print_properties(4)
