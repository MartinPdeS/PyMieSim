"""
Scatterer Footprint Calculation and Visualization
=================================================

This example demonstrates how to compute and visualize the footprint of a scatterer using PyMieSim.
"""

# Import necessary components from PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single.source import Gaussian, PolarizationState
from PyMieSim.single.representations import Footprint

polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=1 * ureg.micrometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=2 * ureg.micrometer,
    source=source,
    medium_refractive_index=1.0 * ureg.RIU,
    refractive_index=1.8 * ureg.RIU,
)

detector = CoherentMode(
    mode_number="HG02",
    numerical_aperture=0.3 * ureg.AU,
    sampling=200,
    gamma_offset=0 * ureg.degree,
    phi_offset=0 * ureg.degree,
    rotation=0 * ureg.degree,
    medium_refractive_index=1.0 * ureg.RIU,
)

footprint = Footprint(
    scatterer=scatterer,
    detector=detector,
    sampling=100,
)

figure = footprint.plot()
