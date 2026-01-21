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

source = Gaussian(
    wavelength=750 * ureg.nanometer,  # 750 nm
    polarization=30 * ureg.degree,  # Right circular polarization
    optical_power=1 * ureg.watt,  # Power in watt
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

scatterer = Cylinder(
    diameter=300 * ureg.nanometer,  # 300 nm
    source=source,
    refractive_index=(1.4 + 0.1j) * ureg.RIU,
    medium_refractive_index=1.33 * ureg.RIU,
)

scatterer.print_properties(4)
