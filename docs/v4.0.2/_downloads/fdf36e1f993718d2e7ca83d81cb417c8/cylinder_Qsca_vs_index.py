"""
InfiniteCylinder: Qsca vs Index
===============================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of refractive index for cylindrical scatterers using PyMieSim, considering multiple wavelengths.
"""
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import InfiniteCylinderSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup

polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[500, 1000, 1500] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = InfiniteCylinderSet(
    diameter=[800] * ureg.nanometer,
    material=np.linspace(1.3, 1.9, 1500) * ureg.RIU,
    medium=[1.0] * ureg.RIU,
)

experiment = Setup(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca", "Qext")

dataframe.plot(x="scatterer:material")
