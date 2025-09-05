"""
Integrating sphere
==================

This example demonstrates the initialization and visualization of an Integrating Sphere detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from TypedUnit import ureg

from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single import plot_system

# %%
# Initializing the detector
detector = IntegratingSphere(
    sampling=500 * ureg.AU,  # Number of sampling points
)

# %%
# Plotting the detector
plot_system(detector)
