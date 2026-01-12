"""
Stokes Parameters Computation
=============================

This example demonstrates the computation and visualization of the Stokes parameters using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

source = Gaussian(
    wavelength=750 * ureg.nanometer,
    polarization=10 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=300 * ureg.nanometer,
    source=source,
    medium_property=1.0 * ureg.RIU,
    property=1.4 * ureg.RIU,
)

data = scatterer.get_stokes(sampling=100)

figure = data.plot()
