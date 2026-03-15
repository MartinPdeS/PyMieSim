"""
InfiniteCylinder: Qsca vs Wavelength
====================================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of wavelength for cylindrical scatterers using PyMieSim, considering cylinders with different diameters and refractive indices.
"""
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import InfiniteCylinderSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup

polarization_state = PolarizationSet(
    angles=0 * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(400, 1000, 150)
    * ureg.nanometer,
    polarization=polarization_state,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = InfiniteCylinderSet(
    diameter=[200, 150] * ureg.nanometer,
    material=[2, 3, 4] * ureg.RIU,
    medium=[1] * ureg.RIU,
)

experiment = Setup(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="source:wavelength", std="scatterer:material")
