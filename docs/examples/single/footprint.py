"""
Scatterer Footprint Calculation and Visualization
=================================================

This example demonstrates how to compute and visualize the footprint of a scatterer using PyMieSim.
"""

# Import necessary components from PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single import Setup

polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=1 * ureg.micrometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=2 * ureg.micrometer,
    medium=1.0 * ureg.RIU,
    material=1.8 * ureg.RIU,
)

detector = CoherentMode(
    mode_number="HG02",
    numerical_aperture=0.3 * ureg.AU,
    sampling=200,
    gamma_offset=0 * ureg.degree,
    phi_offset=0 * ureg.degree,
    rotation=0 * ureg.degree,
    medium=1.0 * ureg.RIU,
)

setup = Setup(
    scatterer=scatterer,
    source=source,
    detector=detector,
)

footprint = setup.get_representation("footprint", sampling=100)

figure = footprint.plot()
