"""
Photodiode Detector
===================

This example demonstrates the initialization and visualization of a Photodiode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import Photodiode
from PyMieSim.units import AU, degree
from PyMieSim.single import plot_system

# %%
# Initializing the detector
detector = Photodiode(
    NA=0.3 * AU,  # Numerical Aperture
    cache_NA=0.2 * AU,
    sampling=500 * AU,  # Number of sampling points
    gamma_offset=0 * degree,  # Gamma offset in degrees
    phi_offset=0 * degree,  # Phi offset in degrees
    polarization_filter=None  # No polarization filter applied
)

# %%
# Plotting the detector
plot_system(detector)
