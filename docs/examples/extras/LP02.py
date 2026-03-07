"""
LP02 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP02 Mode detector using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="LP11",
    sampling=500 * ureg.AU,
    numerical_aperture=0.3 * ureg.AU,
    gamma_offset=0 * ureg.degree,
    rotation=0 * ureg.degree,
    phi_offset=40 * ureg.degree,
)

plotter = SystemPlotter()
plotter.plot(detector)
