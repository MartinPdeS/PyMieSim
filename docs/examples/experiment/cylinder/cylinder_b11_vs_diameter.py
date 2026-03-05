"""
InfiniteCylinder: B1 Scattering Coefficient
===========================================

This example demonstrates how to compute and visualize the B1 scattering coefficient as a function of diameter for cylindrical scatterers using PyMieSim.
"""

import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import InfiniteCylinderSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup

polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(100, 10000, 800) * ureg.nanometer,
    refractive_index=[1.4] * ureg.RIU,
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("b11")

dataframe.plot(x="scatterer:diameter")
