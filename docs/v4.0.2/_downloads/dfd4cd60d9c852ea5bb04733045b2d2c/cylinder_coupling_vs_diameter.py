"""
InfiniteCylinder: Coupling vs Diameter
======================================

This example demonstrates how to compute and visualize the coupling efficiency as a function of diameter for cylindrical scatterers using PyMieSim.
"""

import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import PhotodiodeSet
from PyMieSim.experiment.scatterer import InfiniteCylinderSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup

polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[100, 1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(100, 300, 200) * ureg.nanometer,
    refractive_index=[1.4] * ureg.RIU,
    medium_refractive_index=[1.0] * ureg.RIU,
    source=source,
)

detector = PhotodiodeSet(
    numerical_aperture=[0.1] * ureg.AU,
    phi_offset=[-180.0] * ureg.degree,
    gamma_offset=[0.0] * ureg.degree,
    sampling=[600],
    polarization_filter=None,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="scatterer:diameter")
