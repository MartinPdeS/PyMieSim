"""
CoreShell: Qback vs Core Diameter
=======================================

This example demonstrates how to compute and visualize the backscattering efficiency (Qback)
as functions of core diameter for CoreShell scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import CoreShellSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import print_available, SellmeierMaterial, TabulatedMaterial, SellmeierMedium

print_available()

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[800, 900, 1000] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = CoreShellSet(
    core_diameter=numpy.geomspace(100, 900, 400) * ureg.nanometer,
    shell_thickness=[800] * ureg.nanometer,
    core_material=[TabulatedMaterial("silver")],
    shell_material=[SellmeierMaterial("BK7")],
    medium=[SellmeierMedium("water")],
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qback")

dataframe.plot(x="scatterer:core_diameter")
