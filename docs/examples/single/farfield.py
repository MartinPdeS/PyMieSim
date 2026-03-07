"""
Far-Fields Computation and Visualization
========================================

This example demonstrates the process of computing and visualizing the far-fields of a scatterer using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian, PolarizationState
from PyMieSim.single.representations import FarField

polarization = PolarizationState(angle=30 * ureg.degree)

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=1500 * ureg.nanometer,
    source=source,
    material=1.4 * ureg.RIU,
    medium=1.0 * ureg.RIU,
)

far_fields = FarField(scatterer=scatterer, sampling=100)


figure = far_fields.plot()
