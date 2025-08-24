"""
S1 S2 Function Computation
==========================

This example demonstrates how to compute and visualize the S1 and S2 scattering functions using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from TypedUnit import ureg
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=450 * ureg.nanometer,  # 450 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1 * ureg.watt,  # Arbitrary units
    NA=0.3 * ureg.AU  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=6 * ureg.nanometer,  # 6 nm
    source=source,
    medium_property=1.0 * ureg.RIU,  # Refractive property of the surrounding medium
    property=1.4 * ureg.RIU  # Refractive property of the scatterer
)

# %%
# Computing the data
data = scatterer.get_s1s2(sampling=200)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()
