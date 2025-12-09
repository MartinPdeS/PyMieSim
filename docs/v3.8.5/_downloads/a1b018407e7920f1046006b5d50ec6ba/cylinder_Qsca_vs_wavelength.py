"""
Cylinder: Qsca vs Wavelength
============================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of wavelength for cylindrical scatterers using PyMieSim, considering cylinders with different diameters and refractive indices.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

source = Gaussian(
    wavelength=np.linspace(400, 1000, 150)
    * ureg.nanometer,  # Wavelengths ranging from 400 nm to 1000 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU,  # Numerical Aperture
)

scatterer = Cylinder(
    diameter=[200, 150] * ureg.nanometer,  # Array of diameters: 200 nm, 150 nm, 100 nm
    property=[2, 3, 4] * ureg.RIU,  # Array of refractive indices: 2, 3, 4
    medium_property=[1] * ureg.RIU,  # Refractive index of the surrounding medium
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="source:wavelength", std="scatterer:property")
