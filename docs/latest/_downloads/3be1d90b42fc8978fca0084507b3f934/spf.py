"""
SPF Computation
===============

This example demonstrates the computation and visualization of the Scattering Phase Function (SPF) using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import CoreShell
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
scatterer = CoreShell(
    core_diameter=500 * nanometer,  # 500 nm
    shell_thickness=100 * nanometer,  # 100 nm
    source=source,
    core_property=1.4 * RIU,  # Refractive property of the core
    shell_property=1.8 * RIU,  # Refractive property of the shell
    medium_property=1.0 * RIU,  # Refractive property of the surrounding medium
)

# %%
# Computing the data
data = scatterer.get_spf(sampling=300)  # Specify the number of sampling points

# %%
# Plotting the data
figure = data.plot()
