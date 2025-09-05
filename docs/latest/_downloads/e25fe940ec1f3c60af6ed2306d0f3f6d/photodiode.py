"""
Photodiode Detector
===================

This example demonstrates the initialization and visualization of a Photodiode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from TypedUnit import ureg

from PyMieSim.single.detector import Photodiode
from PyMieSim.single import plot_system

# %%
# Initializing the detector
detector = Photodiode(
    NA=0.3 * ureg.AU,  # Numerical Aperture
    cache_NA=0.2 * ureg.AU,
    sampling=500 * ureg.AU,  # Number of sampling points
    gamma_offset=0 * ureg.degree,  # Gamma offset in degrees
    phi_offset=0 * ureg.degree,  # Phi offset in degrees
)

# %%
# Plotting the detector
plot_system(detector)
