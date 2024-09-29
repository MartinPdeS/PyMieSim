"""
Plot system
===========

Example Script: Using the `plot_system` Function

This script demonstrates how to use the `plot_system` function to create a 3D visualization
of a system consisting of a light source, a scatterer, and a detector. The function leverages
PyVista for rendering the 3D scene.

This script is intended to be used in conjunction with the Read the Docs documentation.
"""

from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import plot_system
from math import sqrt

from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Create a Gaussian light source
source = Gaussian(
    wavelength=1550 * nanometer,  # Wavelength of 1550 nm
    polarization=0 * degree,      # Linear polarization at 0 radians
    optical_power=1 * watt,     # Optical power in arbitrary units
    NA=0.3 * AU               # Numerical Aperture
)

# %%
# Create a Cylinder scatterer
scatterer = Cylinder(
    diameter=7800 * nanometer,     # Diameter of 7.8 micrometers
    source=source,       # The Gaussian source defined above
    medium_index=1.0 * RIU,    # Refractive index of the surrounding medium
    index=sqrt(1.5) * RIU      # Refractive index of the scatterer
)

# %%
# Create a Photodiode detector
detector = CoherentMode(
    mode_number='LP01',
    NA=0.2 * AU,                 # Numerical Aperture
    gamma_offset=0 * degree,        # Gamma offset in degrees
    phi_offset=60 * degree,           # Phi offset in degrees
    polarization_filter=0 * degree   # Polarization filter angle in degrees
)

# %%
# Retrieve the scattering phase function (SPF) from the scatterer
spf = scatterer.get_spf()

# %%
# Visualize the system using the plot_system function
plot_system(source, detector)
