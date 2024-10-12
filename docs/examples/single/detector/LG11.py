"""
Laguerre-Gauss 2-3 Mode Detector
================================

This example demonstrates the initialization and visualization of HG01 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import CoherentMode
from PyMieSim.units import AU, degree

# %%
# Initializing the detector
detector = CoherentMode(
    mode_number="LG23",  # Specifying LP23 mode
    sampling=900 * AU,  # Number of sampling points
    NA=0.4 * AU,  # Numerical Aperture
    gamma_offset=0 * degree,  # Gamma offset
    phi_offset=40 * degree,  # Phi offset in degrees
)

# %%
# Plotting the detector
detector.plot()
