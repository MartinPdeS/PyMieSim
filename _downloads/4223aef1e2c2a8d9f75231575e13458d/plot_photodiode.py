"""
Photodiode Detector
===================

This example demonstrates the initialization and visualization of a Photodiode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import Photodiode

# %%
# Initializing the detector
detector = Photodiode(
    NA=0.3,  # Numerical Aperture
    sampling=500,  # Number of sampling points
    gamma_offset=-45,  # Gamma offset in degrees
    phi_offset=20,  # Phi offset in degrees
    polarization_filter=None  # No polarization filter applied
)

# %%
# Plotting the detector
figure = detector.plot()

# %%
# Displaying the plot
_ = figure.show()
