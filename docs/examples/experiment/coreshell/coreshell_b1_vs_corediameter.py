"""
CoreShell: B1 vs Core Diameter
==============================

This example demonstrates how to compute and visualize the B1 scattering parameter as a function of core diameter for CoreShell scatterers using PyMieSim.
"""

import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import CoreShellSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import Material

polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[800] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = CoreShellSet(
    core_diameter=np.geomspace(100, 3000, 500) * ureg.nanometer,
    shell_thickness=[800] * ureg.nanometer,
    core_refractive_index=[1.3, 1.6] * ureg.RIU,
    shell_material=[Material.BK7],
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("b1")

dataframe.plot(x="scatterer:core_diameter")
