"""
Stokes Parameters Computation
=============================

This example demonstrates the computation and visualization of the Stokes parameters using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from TypedUnit import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
source = Gaussian(
    wavelength=750 * ureg.nanometer,  # 750 nm
    polarization=10 * ureg.degree,  # Right circular polarization
    optical_power=1 * ureg.watt,  # Arbitrary units
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

# %%
# Defining the scatterer
scatterer = Sphere(
    diameter=300 * ureg.nanometer,  # 300 nm
    source=source,
    medium_property=1.0 * ureg.RIU,  # Refractive property of the surrounding medium
    property=1.4 * ureg.RIU,  # Refractive property of the scatterer
)

# %%
# Computing the data
data = scatterer.get_stokes(sampling=100)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()
