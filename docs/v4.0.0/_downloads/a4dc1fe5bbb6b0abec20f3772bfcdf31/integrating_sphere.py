"""
Integrating sphere
==================

This example demonstrates the initialization and visualization of an Integrating Sphere detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single import SystemPlotter

detector = IntegratingSphere(
    sampling=500 * ureg.AU,  # Number of sampling points
)

plotter = SystemPlotter()
plotter.plot(detector)
