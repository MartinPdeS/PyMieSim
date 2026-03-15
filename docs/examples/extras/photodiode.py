"""
Photodiode Detector
===================

This example demonstrates the initialization and visualization of a Photodiode detector using PyMieSim.
"""

import pyvista as pv
from PyMieSim.units import ureg
from PyMieSim.single.detector import Photodiode

detector = Photodiode(
    numerical_aperture=0.3 * ureg.AU,
    cache_numerical_aperture=0.2 * ureg.AU,
    sampling=500 * ureg.AU,
    gamma_offset=45 * ureg.degree,
    phi_offset=0 * ureg.degree,
)

scene = pv.Plotter()

detector.add_to_scene(scene)

scene.show()
