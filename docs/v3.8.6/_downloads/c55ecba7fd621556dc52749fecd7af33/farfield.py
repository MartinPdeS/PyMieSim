"""
Far-Fields Computation and Visualization
========================================

This example demonstrates the process of computing and visualizing the far-fields of a scatterer using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from TypedUnit import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=1000 * ureg.nanometer,  # 1000 nm
    polarization=30 * ureg.degree,  # Right circular polarization
    optical_power=1 * ureg.watt,  # Arbitrary units
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=1500 * ureg.nanometer,  # 1500 nm
    source=source,
    property=1.4 * ureg.RIU,  # Refractive index of the scatterer
    medium_property=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
)

# %%
# Computing the data
data = scatterer.get_farfield(sampling=100)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()
