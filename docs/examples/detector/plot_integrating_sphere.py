"""
Photodiode Detector
===================

This example demonstrates the initialization and visualization of an Integrating Sphere detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import IntegratingSphere

# %%
# Initializing the detector
detector = IntegratingSphere(
    sampling=500,  # Number of sampling points
    polarization_filter=None  # No polarization filter applied
)

# %%
# Plotting the detector
figure = detector.plot()

# %%
# Displaying the plot
_ = figure.show()
