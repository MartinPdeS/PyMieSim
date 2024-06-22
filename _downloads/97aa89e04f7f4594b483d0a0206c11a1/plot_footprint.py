"""
Scatterer Footprint Calculation and Visualization
=================================================

This example demonstrates how to compute and visualize the footprint of a scatterer using PyMieSim.
"""

# Import necessary components from PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single.source import Gaussian
from PyOptik import UsualMaterial

# Define the Gaussian light source with specified properties
source = Gaussian(
    wavelength=1e-6,  # 1000 nm
    polarization=0,
    optical_power=1,  # Arbitrary units
    NA=0.3  # Numerical Aperture
)

# Create a spherical scatterer with a specified diameter and material
scatterer = Sphere(
    diameter=2e-6,  # 2000 nm
    source=source,
    medium_index=1.0,  # Refractive index of the surrounding medium
    material=UsualMaterial.BK7  # Using BK7 glass material
)

# Define the LPMode detector with specific parameters
detector = CoherentMode(
    mode_number="HG02",
    NA=0.3,
    sampling=200,  # Number of sampling points
    gamma_offset=0,
    phi_offset=0,
)

# Compute the footprint data using the defined scatterer and detector
data = detector.get_footprint(scatterer)

# Plot the computed footprint data
figure = data.plot()

# %%
# Display the plot
_ = figure.show()
