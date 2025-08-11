"""
Near-Fields Computation and Visualization
========================================

This example demonstrates the process of computing and visualizing the far-fields of a scatterer using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source
source = Gaussian(
    wavelength=1400 * nanometer,  # 1000 nm
    polarization=30 * degree,  # Right circular polarization
    optical_power=1 * watt,  # Arbitrary units
    NA=0.3 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=1000 * nanometer,  # 1000 nm
    source=source,
    property=1.8 * RIU,  # Refractive index of the scatterer
    medium_property=1.0 * RIU  # Refractive index of the surrounding medium
)

# %%
# Computing the data
data = scatterer.get_nearfield(
    sampling=200,
    field_components=['Ex']
)

# %%
# Plotting the data
figure = data.plot()
