"""
Cylinder: Qsca vs Diameter
==========================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim, considering multiple wavelengths.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

source = Gaussian(
    wavelength=[500, 1000, 1500]
    * ureg.nanometer,  # Array of wavelengths: 500 nm, 1000 nm, 1500 nm
    polarization=30 * ureg.degree,  # Polarization angle in ureg.degrees
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU,  # Numerical Aperture
)

scatterer = Cylinder(
    diameter=np.geomspace(6.36, 10000, 1000)
    * ureg.nanometer,  # Diameters ranging from ~6.36 nm to 10000 nm
    property=[1.4] * ureg.RIU,  # Refractive index of the cylinder
    medium_property=1 * ureg.RIU,  # Refractive index of the surrounding medium
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="scatterer:diameter")
