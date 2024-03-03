"""
Far-Fields Computation and Visualization
========================================

This example demonstrates the process of computing and visualizing the far-fields of a scatterer using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=1e-6,  # 1000 nm
    polarization_value='right',  # Right circular polarization
    polarization_type='circular',
    optical_power=1,  # Arbitrary units
    NA=0.3  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=1.5e-6,  # 1500 nm
    source=source,
    index=1.4,  # Refractive index of the scatterer
    n_medium=1.0  # Refractive index of the surrounding medium
)

# %%
# Computing the data
data = scatterer.get_far_field(sampling=100)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()

# %%
# Display the plot
_ = figure.show()
