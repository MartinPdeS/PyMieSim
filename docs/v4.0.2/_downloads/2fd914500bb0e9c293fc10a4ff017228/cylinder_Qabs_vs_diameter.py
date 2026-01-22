"""
InfiniteCylinder: Qabs vs Diameter
==================================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import InfiniteCylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

source = Gaussian(
    wavelength=400 * ureg.nanometer,  # 400 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU,  # Numerical Aperture
)

scatterer = InfiniteCylinder(
    diameter=np.linspace(1, 800, 300) * ureg.nanometer,  # Diameters ranging from 1 nm to 800 nm
    refractive_index=[
        Material.silver,
        Material.gold,
        Material.aluminium,
    ],  # Scatterer materials
    medium_refractive_index=1 * ureg.RIU,  # Refractive index of the surrounding medium
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get(
    "Qabs"
)  # Assuming Qabs was intended, replace with measure.Qsca if needed

dataframe.plot(x="scatterer:diameter")
