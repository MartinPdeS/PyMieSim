"""
LP11 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP11 Mode detector using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="LP11",
    sampling=300 * ureg.AU,
    numerical_aperture=0.3 * ureg.AU,
    gamma_offset=0 * ureg.degree,
    rotation=0 * ureg.degree,
    phi_offset=30 * ureg.degree,
)

plotter = SystemPlotter()
plotter.plot(detector)
