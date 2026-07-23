"""
InfiniteCylinder: Qabs vs Diameter
==================================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim.
"""
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import InfiniteCylinderSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import print_available, TabulatedMaterial

print_available()


polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(1, 800, 300) * ureg.nanometer,
    material=[TabulatedMaterial("silver"), TabulatedMaterial("gold"), TabulatedMaterial("aluminium")],
    medium=[1],
)

experiment = Setup(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qabs", "Qsca", "Qext")

dataframe.plot(x="scatterer:diameter")
