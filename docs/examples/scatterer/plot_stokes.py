"""
Stokes Parameters Computation
=============================

This example demonstrates the computation and visualization of the Stokes parameters using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=750e-9,  # 750 nm
    polarization_value='right',  # Right circular polarization
    polarization_type='circular',
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

# %%
# Computing the data
data = scatterer.get_stokes(sampling=100)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()

# %%
# Display the plot
_ = figure.show()
