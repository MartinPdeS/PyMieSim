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
    wavelength=1400 * ureg.nanometer,
    polarization=30 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=500 * ureg.nanometer,
    source=source,
    refractive_index=1.8 * ureg.RIU,
    medium_refractive_index=1.0 * ureg.RIU,
)

near_field = NearField(
    scatterer=scatterer,
    sampling=200,
    x_range=(-2 * ureg.micrometer, 2 * ureg.micrometer),
    y_range=(-2 * ureg.micrometer, 2 * ureg.micrometer),
    z=0 * ureg.micrometer,
    field_components=["Hz"]
)

figure = near_field.plot()