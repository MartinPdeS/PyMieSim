"""
Source Plottings
================

This example demonstrates how to visualize the properties of a light source using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=1e-6,  # 1000 nm
    polarization_value=0,  # Linear polarization angle in radians
    polarization_type='linear',
    optical_power=1,  # Arbitrary units
    NA=0.3  # Numerical Aperture
)

# %%
# Plotting the source
figure = source.plot()

# %%
# Display the plot
_ = figure.show()
