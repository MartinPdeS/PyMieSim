"""
SPF Computation
===============

This example demonstrates the computation and visualization of the Scattering Phase Function (SPF) using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=1000e-9,  # 1000 nm
    polarization=0,  # Linear polarization angle in radians
    optical_power=1,  # Arbitrary units
    NA=0.3  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=1.2e-6,  # 1200 nm
    source=source,
    index=1.4,  # Refractive index of the scatterer
    medium_index=1.0,  # Refractive index of the surrounding medium
)

# %%
# Computing the data
data = scatterer.get_spf(sampling=300)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()

# %%
# Display the plot
_ = figure.show()
