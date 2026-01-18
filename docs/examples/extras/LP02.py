"""
LP02 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP02 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="LP11",  # Specifying LP02 mode
    sampling=500 * ureg.AU,  # Number of sampling points
    NA=0.3 * ureg.AU,  # Numerical Aperture
    gamma_offset=0 * ureg.degree,  # Gamma offset
    phi_offset=40 * ureg.degree,  # Phi offset in degrees
)

plotter = SystemPlotter()
plotter.plot(detector)
