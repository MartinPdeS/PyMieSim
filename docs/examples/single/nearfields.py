"""
Near-Fields Computation and Visualization
=========================================

This example demonstrates the process of computing and visualizing the far-fields of a scatterer using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.setup import Setup

polarization_state = PolarizationState(angle=0 * ureg.degree)


source = Gaussian(
    wavelength=300 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=400 * ureg.nanometer,
    material=(1.4 + 5.j) * ureg.RIU,
    medium=1. * ureg.RIU,
)

setup = Setup(
    scatterer=scatterer,
    source=source,
)

near_field = setup.get_representation("nearfields")

near_field.plot(
    "Ex:real",
    "Ex:abs",
    type="scattered",
    plane_origin=(0.0, 0.0, 0.0),
    plane_normal=(0.0, 1.0, 0.0),
    sampling=100,
    extent_scale=4,
    tight_layout=True,
)