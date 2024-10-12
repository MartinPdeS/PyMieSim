"""
S1 S2 Function Computation
==========================

This example demonstrates how to compute and visualize the S1 and S2 scattering functions using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source
source = Gaussian(
    wavelength=450 * nanometer,  # 450 nm
    polarization=0 * degree,  # Linear polarization angle in radians
    optical_power=1 * watt,  # Arbitrary units
    NA=0.3 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=6 * nanometer,  # 6 nm
    source=source,
    medium_property=1.0 * RIU,  # Refractive property of the surrounding medium
    property=1.4 * RIU  # Refractive property of the scatterer
)

# %%
# Computing the data
data = scatterer.get_spf(sampling=200)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()