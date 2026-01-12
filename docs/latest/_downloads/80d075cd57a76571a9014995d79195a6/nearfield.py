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

source = Gaussian(
    wavelength=1400 * ureg.nanometer,
    polarization=30 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=500 * ureg.nanometer,
    source=source,
    property=1.8 * ureg.RIU,
    medium_property=1.0 * ureg.RIU,
)

data = scatterer.get_nearfield(sampling=200, field_components=["Hz"])

figure = data.plot()
