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
from PyMieSim.single.representations import S1S2

source = Gaussian(
    wavelength=450 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=6 * ureg.nanometer,
    source=source,
    medium_refractive_index=1.0 * ureg.RIU,
    refractive_index=1.4 * ureg.RIU,
)

s1s2 = S1S2(scatterer=scatterer, sampling=200)  # Specify the number of sampling points

figure = s1s2.plot()
