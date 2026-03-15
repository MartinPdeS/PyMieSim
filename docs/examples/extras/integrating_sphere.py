"""
Integrating sphere
==================

This example demonstrates the initialization and visualization of an Integrating Sphere detector using PyMieSim.
"""
import pyvista as pv
from PyMieSim.single.detector import IntegratingSphere


detector = IntegratingSphere(
    sampling=500,
)

scene = pv.Plotter()

detector.add_to_scene(scene)

scene.show()
