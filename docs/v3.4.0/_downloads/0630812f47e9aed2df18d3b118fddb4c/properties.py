"""
Print properties
================

This example demonstrates the computation of scattering properties using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyOptik import Material

# %%
# Defining the source
source = Gaussian(
    wavelength=750 * nanometer,  # 750 nm
    polarization=30 * degree,  # Right circular polarization
    optical_power=1 * watt,  # Power in watt
    NA=0.3 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Cylinder(
    diameter=300 * nanometer,  # 300 nm
    source=source,
    property=(1.4 + 0.1j) * RIU,
    medium_property=Material.water
)

scatterer.print_properties()
