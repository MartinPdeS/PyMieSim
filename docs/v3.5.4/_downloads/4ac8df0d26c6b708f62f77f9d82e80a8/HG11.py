"""
Hermite-Gauss 31 Mode Detector
==============================

This example demonstrates the initialization and visualization of HG31 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import CoherentMode
from PyMieSim.units import AU, degree
from PyMieSim.single import plot_system

# %%
# Initializing the detector
detector = CoherentMode(
    mode_number="HG31",  # Specifying HG31 mode
    sampling=500 * AU,  # Number of sampling points
    NA=0.5 * AU,  # Numerical Aperture
    gamma_offset=0 * degree,  # Gamma offset
    phi_offset=40 * degree,  # Phi offset in degrees
    rotation=0 * degree,
)

# %%
# Plotting the detector
plot_system(detector)
