"""
LP11 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP11 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="LP11",  # Specifying LP11 mode
    sampling=300 * ureg.AU,  # Number of sampling points
    NA=0.3 * ureg.AU,  # Numerical Aperture
    gamma_offset=0 * ureg.degree,  # Gamma offset
    rotation=0 * ureg.degree,  # Rotation angle
    phi_offset=30 * ureg.degree,  # Phi offset in degrees
)

plotter = SystemPlotter()
plotter.plot(detector)
