"""
InfiniteCylinder: Coupling vs Wavelength
========================================

This example demonstrates how to compute and visualize the coupling efficiency as a function of wavelength for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import CoherentModeSet
from PyMieSim.experiment.scatterer import InfiniteCylinderSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import Material

polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(950, 1050, 300) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(300, 500, 5) * ureg.nanometer,
    material=[Material.BK7],
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

detector = CoherentModeSet(
    mode_number=["LP11"],
    numerical_aperture=[0.05, 0.01] * ureg.AU,
    phi_offset=[-180] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[300],
    rotation=[0] * ureg.degree,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling", scale_unit=True)

dataframe.plot(x="source:wavelength", std="scatterer:diameter")
