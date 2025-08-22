"""
Scatterer Footprint Calculation and Visualization
=================================================

This example demonstrates how to compute and visualize the footprint of a scatterer using PyMieSim.
"""

# Import necessary components from PyMieSim
from TypedUnit import ureg
from PyOptik import Material

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single.source import Gaussian

# Define the Gaussian light source with specified properties
source = Gaussian(
    wavelength=1 * ureg.micrometer,  # 1000 nm
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,  # Arbitrary units
    NA=0.3 * ureg.AU  # Numerical Aperture
)

# Create a spherical scatterer with a specified diameter and property
scatterer = Sphere(
    diameter=2 * ureg.micrometer,  # 2000 nm
    source=source,
    medium_property=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
    property=Material.BK7  # Using BK7 glass property
)

# Define the LPMode detector with specific parameters
detector = CoherentMode(
    mode_number="HG02",
    NA=0.3 * ureg.AU,
    sampling=200 * ureg.AU,  # Number of sampling points
    gamma_offset=0 * ureg.degree,
    phi_offset=0 * ureg.degree,
)

# Compute the footprint data using the defined scatterer and detector
data = detector.get_footprint(scatterer)

# Plot the computed footprint data
figure = data.plot()
