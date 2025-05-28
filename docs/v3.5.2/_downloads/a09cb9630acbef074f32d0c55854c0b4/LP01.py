"""
LP01 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP01 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import CoherentMode
from PyMieSim.units import AU, degree, RIU
from PyMieSim.single import plot_system

# %%
# Initializing the detector
detector = CoherentMode(
    mode_number="LP01",  # Specifying LP01 mode
    sampling=500 * AU,  # Number of sampling points
    NA=1.0 * AU,  # Numerical Aperture
    gamma_offset=90 * degree,  # Gamma offset
    phi_offset=0 * degree,  # Phi offset in degrees
    medium_refractive_index=1.3 * RIU
)

# %%
# Plotting the detector
plot_system(detector)
