"""
Source Plottings
================

This example demonstrates how to visualize the properties of a light source using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.source import Gaussian
from PyMieSim.units import nanometer, degree, watt, AU

# %%
# Defining the source
source = Gaussian(
    wavelength=1 * nanometer,  # 1000 nm
    polarization=0 * degree,  # Linear polarization angle in radians
    optical_power=1 * watt,  # Arbitrary units
    NA=0.3 * AU  # Numerical Aperture
)

# %%
# Plotting the source
source.plot()
