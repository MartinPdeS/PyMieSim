"""
Scatterer Footprint Calculation and Visualization
=================================================

This example demonstrates how to compute and visualize the footprint of a scatterer using PyMieSim.
"""

# Import necessary components from PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single.source import Gaussian
from PyOptik import Material
from PyMieSim.units import micrometer, degree, watt, AU, RIU

# Define the Gaussian light source with specified properties
source = Gaussian(
    wavelength=1 * micrometer,  # 1000 nm
    polarization=0 * degree,
    optical_power=1 * watt,  # Arbitrary units
    NA=0.3 * AU  # Numerical Aperture
)

# Create a spherical scatterer with a specified diameter and property
scatterer = Sphere(
    diameter=2 * micrometer,  # 2000 nm
    source=source,
    medium_property=1.0 * RIU,  # Refractive index of the surrounding medium
    property=Material.BK7  # Using BK7 glass property
)

# Define the LPMode detector with specific parameters
detector = CoherentMode(
    mode_number="HG02",
    NA=0.3 * AU,
    sampling=200 * AU,  # Number of sampling points
    gamma_offset=0 * degree,
    phi_offset=0 * degree,
)

# Compute the footprint data using the defined scatterer and detector
data = detector.get_footprint(scatterer)

# Plot the computed footprint data
figure = data.plot()
