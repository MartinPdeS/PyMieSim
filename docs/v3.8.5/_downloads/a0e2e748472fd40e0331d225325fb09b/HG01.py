"""
Hermite-Gauss 01 Mode Detector
==============================

This example demonstrates the initialization and visualization of HG01 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="HG01",  # Specifying LP01 mode
    sampling=500 * ureg.AU,  # Number of sampling points
    NA=0.5 * ureg.AU,  # Numerical Aperture
    rotation=0 * ureg.degree,  # Rotation angle
    gamma_offset=0 * ureg.degree,  # Gamma offset
    phi_offset=40 * ureg.degree,  # Phi offset in degrees
)

plotter = SystemPlotter()
plotter.plot(detector)
