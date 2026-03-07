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

from PyMieSim.binary.interface_dispersive_material import DispersiveMaterial


material_bk7 = DispersiveMaterial(
    name="BK7",
    wavelengths=[400, 500, 600, 700, 800] * ureg.nanometer,
    refractive_indices=[1.5302, 1.5195, 1.5143, 1.5114, 1.5098] * ureg.RIU,
)

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
    medium=1.0 * ureg.RIU,
    material=material_bk7,
)

scatterer.Qsca

s1s2 = S1S2(scatterer=scatterer, sampling=200)  # Specify the number of sampling points

figure = s1s2.plot()
