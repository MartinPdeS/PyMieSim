"""
Scatterer Footprint Calculation and Visualization
=================================================

This example demonstrates how to compute and visualize the footprint of a scatterer using PyMieSim.
"""

# Import necessary components from PyMieSim
from PyMieSim.units import ureg
from PyOptik import Material

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.detector import CoherentMode
from PyMieSim.single.source import Gaussian

source = Gaussian(
    wavelength=1 * ureg.micrometer,
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=2 * ureg.micrometer,
    source=source,
    medium_property=1.0 * ureg.RIU,
    property=Material.BK7,
)

detector = CoherentMode(
    mode_number="HG02",
    NA=0.3 * ureg.AU,
    sampling=200 * ureg.AU,
    gamma_offset=0 * ureg.degree,
    phi_offset=0 * ureg.degree,
)

data = detector.get_footprint(scatterer)

figure = data.plot()
