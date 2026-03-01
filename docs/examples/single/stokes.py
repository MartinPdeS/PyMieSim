"""
Stokes Parameters Computation
=============================

This example demonstrates the computation and visualization of the Stokes parameters using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian, PolarizationState
from PyMieSim.single.representations import Stokes

polarization_state = PolarizationState(angle=10 * ureg.degree)

source = Gaussian(
    wavelength=750 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=300 * ureg.nanometer,
    source=source,
    medium_refractive_index=1.0 * ureg.RIU,
    refractive_index=1.4 * ureg.RIU,
)

stokes = Stokes(scatterer=scatterer, sampling=100)

figure = stokes.plot()
