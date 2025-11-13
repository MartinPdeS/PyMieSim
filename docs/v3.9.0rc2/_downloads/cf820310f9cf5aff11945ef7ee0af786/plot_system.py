"""
Plot system
===========

Example Script: Using the `plot_system` Function

This script demonstrates how to use the `plot_system` function to create a 3D visualization
of a system consisting of a light source, a scatterer, and a detector. The function leverages
PyVista for rendering the 3D scene.

This script is intended to be used in conjunction with the Read the Docs documentation.
"""

from math import sqrt
from TypedUnit import ureg

from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import plot_system

# %%
# Create a Gaussian light source
source = Gaussian(
    wavelength=1550 * ureg.nanometer,  # Wavelength of 1550 nm
    polarization=0 * ureg.degree,  # Linear polarization at 0 radians
    optical_power=1 * ureg.watt,  # Optical power in arbitrary units
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

# %%
# Create a Cylinder scatterer
scatterer = Cylinder(
    diameter=7800 * ureg.nanometer,  # Diameter of 7.8 micrometers
    source=source,  # The Gaussian source defined above
    medium_property=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
    property=sqrt(1.5) * ureg.RIU,  # Refractive index of the scatterer
)

# %%
# Create a Photodiode detector
detector = CoherentMode(
    mode_number="LP01",
    NA=0.2 * ureg.AU,  # Numerical Aperture
    gamma_offset=0 * ureg.degree,  # Gamma offset in degrees
    phi_offset=60 * ureg.degree,  # Phi offset in degrees
    polarization_filter=0 * ureg.degree,  # Polarization filter angle in degrees
)

# %%
# Retrieve the scattering phase function (SPF) from the scatterer
spf = scatterer.get_spf()

# %%
# Visualize the system using the plot_system function
plot_system(source, detector)
