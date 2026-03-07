"""
LP01 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP01 Mode detector using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter

detector = CoherentMode(
    mode_number="LP01",
    sampling=500 * ureg.AU,
    numerical_aperture=1.0 * ureg.AU,
    gamma_offset=90 * ureg.degree,
    rotation=0 * ureg.degree,
    phi_offset=0 * ureg.degree,
    medium_refractive_index=1.3 * ureg.RIU,
)

plotter = SystemPlotter()
plotter.plot(detector)
