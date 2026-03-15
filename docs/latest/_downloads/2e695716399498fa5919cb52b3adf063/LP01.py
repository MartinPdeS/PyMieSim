"""
LP01 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP01 Mode detector using PyMieSim.
"""

import pyvista as pv
from PyMieSim.units import ureg
from PyMieSim.single.detector import CoherentMode


detector = CoherentMode(
    mode_number="LP01",
    sampling=500 * ureg.AU,
    numerical_aperture=1.0 * ureg.AU,
    gamma_offset=90 * ureg.degree,
    rotation=0 * ureg.degree,
    phi_offset=0 * ureg.degree,
    medium=1.0 * ureg.RIU,
)

scene = pv.Plotter()

detector.add_to_scene(scene)

scene.show()
