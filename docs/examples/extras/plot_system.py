"""
Example Script: Using the `plot_system` Function

This script demonstrates how to use the `plot_system` function to create a 3D visualization
of a system consisting of a light source, a scatterer, and a detector. The function leverages
PyVista for rendering the 3D scene.

This script is intended to be used in conjunction with the Read the Docs documentation.
"""

from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import plot_system
from math import sqrt

# %%
# Create a Gaussian light source
source = Gaussian(
    wavelength=1550e-9,  # Wavelength of 1550 nm
    polarization=0,      # Linear polarization at 0 radians
    optical_power=1,     # Optical power in arbitrary units
    NA=0.3               # Numerical Aperture
)

# %%
# Create a Cylinder scatterer
scatterer = Cylinder(
    diameter=7.8e-6,     # Diameter of 7.8 micrometers
    source=source,       # The Gaussian source defined above
    medium_index=1.0,    # Refractive index of the surrounding medium
    index=sqrt(1.5)      # Refractive index of the scatterer
)

# %%
# Create a Photodiode detector
detector = Photodiode(
    NA=0.1,                 # Numerical Aperture
    gamma_offset=90,        # Gamma offset in degrees
    phi_offset=0,           # Phi offset in degrees
    polarization_filter=0   # Polarization filter angle in degrees
)

# %%
# Retrieve the scattering phase function (SPF) from the scatterer
spf = scatterer.get_spf()

# %%
# Visualize the system using the plot_system function
plot_system(source, detector, spf)
