"""
InfiniteCylinder: Qsca vs Index
===============================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of refractive index for cylindrical scatterers using PyMieSim, considering multiple wavelengths.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import InfiniteCylinderSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup

polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[500, 1000, 1500] * ureg.nanometer,  # Array of wavelengths: 500 nm, 1000 nm, 1500 nm
    polarization=polarization_set,  # Polarization angle in ureg.degrees
    optical_power=[1e-3] * ureg.watt,  # 1 milliureg.watt
    numerical_aperture=[0.2] * ureg.AU,  # Numerical Aperture
)

scatterer = InfiniteCylinderSet(
    diameter=[800] * ureg.nanometer,  # Fixed diameter of 800 nm
    refractive_index=np.linspace(1.3, 1.9, 1500) * ureg.RIU,  # Refractive index ranging from 1.3 to 1.9
    medium_refractive_index=[1.0] * ureg.RIU,  # Refractive index of the surrounding medium
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca", "Qext")

dataframe.plot(x="scatterer:refractive_index")
