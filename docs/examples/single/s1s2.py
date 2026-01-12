"""
S1 S2 Function Computation
==========================

This example demonstrates how to compute and visualize the S1 and S2 scattering functions using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

source = Gaussian(
    wavelength=450 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=6 * ureg.nanometer,
    source=source,
    medium_property=1.0 * ureg.RIU,
    property=1.4 * ureg.RIU,
)

data = scatterer.get_s1s2(sampling=200)  # Specify the number of sampling points

figure = data.plot()
