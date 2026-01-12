"""
Source Plottings
================

This example demonstrates how to visualize the properties of a light source using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.source import Gaussian

source = Gaussian(
    wavelength=1 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

source.plot()
