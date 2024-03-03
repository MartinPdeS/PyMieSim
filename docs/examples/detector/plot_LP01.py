"""
LP01 Mode Detector
==================

This example demonstrates the initialization and visualization of an LP01 Mode detector using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import LPMode

# %%
# Initializing the detector
detector = LPMode(
    mode_number="LP01",  # Specifying LP01 mode
    sampling=500,  # Number of sampling points
    NA=0.5,  # Numerical Aperture
    gamma_offset=0,  # Gamma offset
    phi_offset=40,  # Phi offset in degrees
    coupling_mode='Point'  # Coupling mode set to 'Point'
)

# %%
# Plotting the detector
figure = detector.plot()

# %%
# Customizing the plot appearance
figure.background_color = 'black'  # Setting the background color to black
figure.unit_size = (1200, 1200)  # Adjusting the size of the plot

# %%
# Displaying the plot
_ = figure.show()
