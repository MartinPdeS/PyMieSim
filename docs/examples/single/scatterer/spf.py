"""
SPF Computation
===============

This example demonstrates the computation and visualization of the Scattering Phase Function (SPF) using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source
source = Gaussian(
    wavelength=1000 * nanometer,  # 1000 nm
    polarization=0 * degree,  # Linear polarization angle in radians
    optical_power=1 * watt,  # Arbitrary units
    NA=0.3 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=2000 * nanometer,  # 1200 nm
    source=source,
    property=1.4 * RIU,  # Refractive property of the scatterer
    medium_property=1.0 * RIU,  # Refractive property of the surrounding medium
)

# %%
# Computing the data
data = scatterer.get_spf(sampling=300)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()