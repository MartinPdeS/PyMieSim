"""
Laguerre-Gauss 2-3 Mode Detector
================================

This example demonstrates the initialization and visualization of HG01 Mode detector using PyMieSim.
"""

import pyvista as pv
from PyMieSim.units import ureg
from PyMieSim.single.detector import CoherentMode


detector = CoherentMode(
    mode_number="LG23",
    sampling=900 * ureg.AU,
    numerical_aperture=0.4 * ureg.AU,
    gamma_offset=0 * ureg.degree,
    rotation=0 * ureg.degree,
    phi_offset=40 * ureg.degree,
)

scene = pv.Plotter()

detector.add_to_scene(scene)

scene.show()
