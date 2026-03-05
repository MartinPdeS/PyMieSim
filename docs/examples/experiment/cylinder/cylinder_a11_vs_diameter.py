"""
InfiniteCylinder: A1 Scattering Coefficient
===========================================

This example demonstrates how to compute and visualize the A1 scattering coefficient as a function of diameter for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import InfiniteCylinderSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup

polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,  # 400 nm
    polarization=polarization_set,  # Polarization angle in ureg.degrees
    optical_power=[1e-3] * ureg.watt,  # 1 milliureg.watt
    numerical_aperture=[0.2] * ureg.AU,  # Numerical Aperture
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(100, 10000, 800)
    * ureg.nanometer,  # Diameters ranging from 100 nm to 10000 nm
    refractive_index=[1.4] * ureg.RIU,  # Refractive index of the cylinder
    medium_refractive_index=[1.0] * ureg.RIU,  # Refractive index of the surrounding medium
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("a12")

dataframe.plot(x="scatterer:diameter")
