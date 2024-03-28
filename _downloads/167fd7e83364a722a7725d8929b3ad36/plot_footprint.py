"""
Scatterer Footprint Calculation and Visualization
=================================================

This example demonstrates how to compute and visualize the footprint of a scatterer using PyMieSim.
"""

# Import necessary components from PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.detector import LPMode
from PyMieSim.single.source import Gaussian
from PyMieSim.materials import BK7

# Define the Gaussian light source with specified properties
source = Gaussian(
    wavelength=1e-6,  # 1000 nm
    polarization_value=0,
    polarization_type='linear',
    optical_power=1,  # Arbitrary units
    NA=0.3  # Numerical Aperture
)

# Create a spherical scatterer with a specified diameter and material
scatterer = Sphere(
    diameter=2e-6,  # 2000 nm
    source=source,
    material=BK7  # Using BK7 glass material
)

# Define the LPMode detector with specific parameters
detector = LPMode(
    mode_number="LP21",
    NA=0.3,
    sampling=200,  # Number of sampling points
    gamma_offset=0,
    phi_offset=0,
    coupling_mode='Point'  # Coupling mode
)

# Compute the footprint data using the defined scatterer and detector
data = detector.get_footprint(scatterer)

# Plot the computed footprint data
figure = data.plot()

# %%
# Display the plot
_ = figure.show()
