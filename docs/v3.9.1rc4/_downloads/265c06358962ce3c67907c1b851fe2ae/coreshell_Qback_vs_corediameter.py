"""
CoreShell: Qback vs Core Diameter
=======================================

This example demonstrates how to compute and visualize the backscattering efficiency (Qback)
as functions of core diameter for CoreShell scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

source = Gaussian(
    wavelength=[800, 900, 1000]
    * ureg.nanometer,  # Array of wavelengths: 800 nm, 900 nm, 1000 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU,  # Numerical Aperture
)

scatterer = CoreShell(
    core_diameter=numpy.geomspace(100, 600, 400)
    * ureg.nanometer,  # Core diameters from 100 nm to 600 nm
    shell_thickness=800 * ureg.nanometer,  # Shell width of 800 nm
    core_property=Material.silver,  # Core material
    shell_property=Material.BK7,  # Shell material
    medium_property=1 * ureg.RIU,  # Surrounding medium's refractive index
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qback")

dataframe.plot(x="scatterer:core_diameter")
