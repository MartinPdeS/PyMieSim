"""
InfiniteCylinder: Qabs vs Diameter
==================================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import InfiniteCylinderSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import Material

polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(1, 800, 300) * ureg.nanometer,
    material=[Material.silver, Material.gold, Material.aluminium],
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qabs")

dataframe.plot(x="scatterer:diameter")
