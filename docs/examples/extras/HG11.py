"""
Hermite-Gauss 31 Mode Detector
==============================

This example demonstrates the initialization and visualization of HG31 Mode detector using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="HG31",
    sampling=500 * ureg.AU,
    numerical_aperture=0.5 * ureg.AU,
    gamma_offset=0 * ureg.degree,
    phi_offset=40 * ureg.degree,
    rotation=0 * ureg.degree,
)

plotter = SystemPlotter()
plotter.plot(detector)
