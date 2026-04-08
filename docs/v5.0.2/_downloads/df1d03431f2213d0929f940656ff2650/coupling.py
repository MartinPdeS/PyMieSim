"""
Coupling
========

This example demonstrates the process of computing the coupling efficiency of a scatterer to a detector using PyMieSim.
"""

from PyMieSim.units import ureg
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.detector import Photodiode
from PyMieSim.single.setup import Setup


polarization = PolarizationState(angle=30 * ureg.degree)

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = Sphere(
    diameter=1500 * ureg.nanometer,
    material=1.45,
    medium=1.0,
)

detector = Photodiode(
    numerical_aperture=0.3,
    phi_offset=45 * ureg.degree,
    gamma_offset=45 * ureg.degree,
)


setup = Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

value = setup.get("coupling")

print(value)
