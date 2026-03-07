"""
Stokes Parameters Computation
=============================

This example demonstrates the computation and visualization of the Stokes parameters using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import PolarizationState
from PyMieSim.single.representations import Stokes

polarization_state = PolarizationState(angle=10 * ureg.degree)

source = Gaussian(
    wavelength=750 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=300 * ureg.nanometer,
    source=source,
    medium=1.0 * ureg.RIU,
    material=1.4 * ureg.RIU,
)

stokes = Stokes(scatterer=scatterer, sampling=100)

figure = stokes.plot()
