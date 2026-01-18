"""
Print properties
================

This example demonstrates the computation of scattering properties using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyOptik import Material

source = Gaussian(
    wavelength=750 * ureg.nanometer,  # 750 nm
    polarization=30 * ureg.degree,  # Right circular polarization
    optical_power=1 * ureg.watt,  # Power in watt
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

scatterer = Cylinder(
    diameter=300 * ureg.nanometer,  # 300 nm
    source=source,
    property=(1.4 + 0.1j) * ureg.RIU,
    medium_property=Material.water,
)

scatterer.print_properties()
