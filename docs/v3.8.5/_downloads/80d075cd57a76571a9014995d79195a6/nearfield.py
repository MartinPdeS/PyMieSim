"""
Near-Fields Computation and Visualization
=========================================

This example demonstrates the process of computing and visualizing the far-fields of a scatterer using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from TypedUnit import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

source = Gaussian(
    wavelength=1400 * ureg.nanometer,  # 1000 nm
    polarization=30 * ureg.degree,  # Right circular polarization
    optical_power=1 * ureg.watt,  # Arbitrary units
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

scatterer = Sphere(
    diameter=500 * ureg.nanometer,  # 500 nm
    source=source,
    property=1.8 * ureg.RIU,  # Refractive index of the scatterer
    medium_property=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
)

data = scatterer.get_nearfield(sampling=200, field_components=["Hz"])

figure = data.plot()
