"""
Laguerre-Gauss 2-3 Mode Detector
================================

This example demonstrates the initialization and visualization of HG01 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="LG23",  # Specifying LP23 mode
    sampling=900 * ureg.AU,  # Number of sampling points
    NA=0.4 * ureg.AU,  # Numerical Aperture
    gamma_offset=0 * ureg.degree,  # Gamma offset
    rotation=0 * ureg.degree,  # Rotation angle
    phi_offset=40 * ureg.degree,  # Phi offset in degrees
)

plotter = SystemPlotter()
plotter.plot(detector)
