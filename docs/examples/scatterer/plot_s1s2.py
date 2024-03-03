"""
S1 S2 Function Computation
==========================

This example demonstrates how to compute and visualize the S1 and S2 scattering functions using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=450e-9,  # 450 nm
    polarization_value=0,  # Linear polarization angle in radians
    polarization_type='linear',
    optical_power=1,  # Arbitrary units
    NA=0.3  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=6e-9,  # 6 nm
    source=source,
    index=1.4  # Refractive index of the scatterer
)

# %%
# Computing the data
data = scatterer.get_s1s2(sampling=200)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()

# %%
# Display the plot
_ = figure.show()
