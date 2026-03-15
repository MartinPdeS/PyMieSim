"""
Samples Properties
==================
"""

from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState

polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=800 * ureg.nanometer,
    material=1.4 * ureg.RIU,
    medium=1.0 * ureg.RIU,
)

scatterer.init(source)

scatterer.print_properties(4)
