"""
Scattering properties
=====================

This example demonstrates the computation of scattering properties using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source
source = Gaussian(
    wavelength=750 * nanometer,  # 750 nm
    polarization=30 * degree,  # Right circular polarization
    optical_power=1 * watt,  # Arbitrary units
    NA=0.3 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=300 * nanometer,  # 300 nm
    source=source,
    index=1.4 * RIU,  # Refractive index of the scatterer
    medium_index=1 * RIU
)

scatterer.print_properties()
