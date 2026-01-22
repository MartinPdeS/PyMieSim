"""
Near-Fields Computation and Visualization
=========================================

This example demonstrates the process of computing and visualizing the far-fields of a scatterer using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.representations import NearField

source = Gaussian(
    wavelength=300 * ureg.nanometer,
    polarization=90 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=30 * ureg.nanometer,
    source=source,
    refractive_index=(1.4 + 0.j) * ureg.RIU,
    medium_refractive_index=1. * ureg.RIU,
)

near_field = NearField(
    scatterer=scatterer,
)

near_field.plot(
    "|E|:abs",
    type="scattered",
    plane_origin=(0.0, 0.0, 0.0),
    plane_normal=(0.0, 1.0, 0.0),
    sampling=400,
    extent_scale=100,
)