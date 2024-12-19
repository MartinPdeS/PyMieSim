"""
Integrating sphere
==================

This example demonstrates the initialization and visualization of an Integrating Sphere detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.units import AU

# %%
# Initializing the detector
detector = IntegratingSphere(
    sampling=500 * AU,  # Number of sampling points
    polarization_filter=None  # No polarization filter applied
)

# %%
# Plotting the detector
detector.plot()