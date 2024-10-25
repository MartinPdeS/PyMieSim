"""
Photodiode Detector
===================

This example demonstrates the initialization and visualization of a Photodiode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import Photodiode
from PyMieSim.units import AU, degree

# %%
# Initializing the detector
detector = Photodiode(
    NA=0.2 * AU,  # Numerical Aperture
    cache_NA=0.1 * AU,  # Numerical Aperture
    sampling=400 * AU,  # Number of sampling points
    gamma_offset=0 * degree,  # Gamma offset in degrees
    phi_offset=0 * degree,  # Phi offset in degrees
    polarization_filter=None  # No polarization filter applied
)

# %%
# Plotting the detector
detector.plot()
