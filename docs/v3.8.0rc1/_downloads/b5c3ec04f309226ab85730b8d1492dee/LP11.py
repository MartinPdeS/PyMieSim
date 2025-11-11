"""
LP11 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP11 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from TypedUnit import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import plot_system

# %%
# Initializing the detector
detector = CoherentMode(
    mode_number="LP11",  # Specifying LP11 mode
    sampling=300 * ureg.AU,  # Number of sampling points
    NA=0.3 * ureg.AU,  # Numerical Aperture
    gamma_offset=0 * ureg.degree,  # Gamma offset
    phi_offset=30 * ureg.degree,  # Phi offset in degrees
)

# %%
# Plotting the detector
plot_system(detector)
