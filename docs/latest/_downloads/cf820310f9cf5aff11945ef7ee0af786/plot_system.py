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
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter
from PyMieSim.single.representations import SPF

source = Gaussian(
    wavelength=1550 * ureg.nanometer,  # Wavelength of 1550 nm
    polarization=0 * ureg.degree,  # Linear polarization at 0 radians
    optical_power=1 * ureg.watt,  # Optical power in arbitrary units
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

scatterer = Cylinder(
    diameter=7800 * ureg.nanometer,  # Diameter of 7.8 micrometers
    source=source,  # The Gaussian source defined above
    medium_refractive_index=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
    refractive_index=sqrt(1.5) * ureg.RIU,  # Refractive index of the scatterer
)

detector = CoherentMode(
    mode_number="LP01",
    NA=0.2 * ureg.AU,  # Numerical Aperture
    gamma_offset=0 * ureg.degree,  # Gamma offset in degrees
    phi_offset=60 * ureg.degree,  # Phi offset in degrees
    rotation=0 * ureg.degree,  # Rotation angle in degrees
    polarization_filter=0 * ureg.degree,  # Polarization filter angle in degrees
)

spf = SPF(scatterer=scatterer)

plotter = SystemPlotter(show_axis_label=False)

plotter.plot(source, scatterer, detector, data=spf)

