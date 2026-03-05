"""
S1 S2 Function Computation
==========================

This example demonstrates how to compute and visualize the S1 and S2 scattering functions using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian, PolarizationState
from PyMieSim.single.representations import S1S2

polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=450 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=6 * ureg.nanometer,
    source=source,
    medium_refractive_index=1.0 * ureg.RIU,
    refractive_index=1.4 * ureg.RIU,
)

scatterer.Qsca

s1s2 = S1S2(scatterer=scatterer, sampling=200)  # Specify the number of sampling points

figure = s1s2.plot()
