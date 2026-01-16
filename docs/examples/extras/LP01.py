"""
LP01 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP01 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="LP01",  # Specifying LP01 mode
    sampling=500 * ureg.AU,  # Number of sampling points
    NA=1.0 * ureg.AU,  # Numerical Aperture
    gamma_offset=90 * ureg.degree,  # Gamma offset
    phi_offset=0 * ureg.degree,  # Phi offset in degrees
    medium_refractive_index=1.3 * ureg.RIU,
)

plotter = SystemPlotter()
plotter.plot(detector)
