"""
S1 S2 Function Computation
==========================

This example demonstrates how to compute and visualize the S1 and S2 scattering functions using PyMieSim.
"""

from PyMieSim.units import ureg
from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import PolarizationState
from PyMieSim.single.material import Material
from PyMieSim.single import Setup

bk7 = Material(
    name="BK7",
    refractive_indices=[1.3, 1.4, 1.5, 1.6, 1.7] * ureg.RIU,
    wavelengths=[400, 500, 600, 700, 800] * ureg.nanometer,
)


polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=450 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = InfiniteCylinder(
    diameter=6 * ureg.nanometer,
    medium=1.0 * ureg.RIU,
    # material=bk7,
    material=1.4 * ureg.RIU,
)

setup = Setup(
    scatterer=scatterer,
    source=source,
)

s1s2 = setup.get_representation("s1s2", sampling=200)  # Specify the number of sampling points


# s1s2 = S1S2(setup=setup, sampling=200)  # Specify the number of sampling points

figure = s1s2.plot()
