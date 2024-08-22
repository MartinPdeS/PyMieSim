"""
Scattering properties
=====================

This example demonstrates the computation of scattering properties using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=750e-9,  # 750 nm
    polarization=30,  # Right circular polarization
    optical_power=1,  # Arbitrary units
    NA=0.3  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=300e-9,  # 300 nm
    source=source,
    index=1.4  # Refractive index of the scatterer
)

scatterer.print_properties()
