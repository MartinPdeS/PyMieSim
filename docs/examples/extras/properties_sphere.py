"""
Samples Properties
==================
"""

from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.representations import FarField

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=800 * ureg.nanometer,
    source=source,
    refractive_index=1.4 * ureg.RIU,
    medium_refractive_index=1.0 * ureg.RIU,
)

scatterer.print_properties(4)
