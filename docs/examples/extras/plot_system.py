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

from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import PolarizationState
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single import SystemPlotter
from PyMieSim.single.representations import SPF

polarization_state = PolarizationState(
    angle=0 * ureg.degree,  # Linear polarization at 0 degrees
)

source = Gaussian(
    wavelength=1550 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = InfiniteCylinder(
    diameter=7800 * ureg.nanometer,
    source=source,
    medium=1.0 * ureg.RIU,
    material=sqrt(1.5) * ureg.RIU,
)

detector = CoherentMode(
    mode_number="LP01",
    numerical_aperture=0.2 * ureg.AU,
    gamma_offset=0 * ureg.degree,
    phi_offset=60 * ureg.degree,
    rotation=0 * ureg.degree,
    polarization_filter=0 * ureg.degree,
)

spf = SPF(scatterer=scatterer)

plotter = SystemPlotter(show_axis_label=False)

plotter.plot(source, scatterer, detector, data=spf)

