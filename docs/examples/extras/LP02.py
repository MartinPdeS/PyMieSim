"""
LP02 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP02 Mode detector using PyMieSim.
"""

import pyvista as pv
from PyMieSim.units import ureg
from PyMieSim.single.detector import CoherentMode


detector = CoherentMode(
    mode_number="LP11",
    sampling=500,
    numerical_aperture=0.3,
    gamma_offset=0 * ureg.degree,
    rotation=0 * ureg.degree,
    phi_offset=40 * ureg.degree,
)

scene = pv.Plotter()

detector.add_to_scene(scene)

scene.show()