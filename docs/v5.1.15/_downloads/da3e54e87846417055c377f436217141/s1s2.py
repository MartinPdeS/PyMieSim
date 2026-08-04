"""
S1 S2 Function Computation
==========================

This example demonstrates how to compute and visualize the S1 and S2 scattering functions using PyMieSim.
"""

from PyMieSim.units import ureg
from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.material import SellmeierMaterial
from PyMieSim.single.setup import Setup

bk7 = SellmeierMaterial("BK7")

polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=450 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = InfiniteCylinder(
    diameter=6 * ureg.nanometer,
    medium=1.0,
    material=bk7,
)

setup = Setup(
    scatterer=scatterer,
    source=source,
)

s1s2 = setup.get_representation("s1s2", sampling=200)

figure = s1s2.plot()
