"""
Samples Properties
==================
"""

from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian, PolarizationState

polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=800 * ureg.nanometer,
    source=source,
    refractive_index=1.4 * ureg.RIU,
    medium_refractive_index=1.0 * ureg.RIU,
)

scatterer.print_properties(4)
