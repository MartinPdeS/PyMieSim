"""
Hermite-Gauss 31 Mode Detector
==============================

This example demonstrates the initialization and visualization of HG31 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import CoherentMode

# %%
# Initializing the detector
detector = CoherentMode(
    mode_number="HG31",  # Specifying HG31 mode
    sampling=500,  # Number of sampling points
    NA=0.5,  # Numerical Aperture
    gamma_offset=0,  # Gamma offset
    phi_offset=40,  # Phi offset in degrees
    rotation=0,
)

# %%
# Plotting the detector
detector.plot()
