"""
Laguerre-Gauss 2-3 Mode Detector
================================

This example demonstrates the initialization and visualization of HG01 Mode detector using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="LG23",
    sampling=900 * ureg.AU,
    numerical_aperture=0.4 * ureg.AU,
    gamma_offset=0 * ureg.degree,
    rotation=0 * ureg.degree,
    phi_offset=40 * ureg.degree,
)

plotter = SystemPlotter()
plotter.plot(detector)
