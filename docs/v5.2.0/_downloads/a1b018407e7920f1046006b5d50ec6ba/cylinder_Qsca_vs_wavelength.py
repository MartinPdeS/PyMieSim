"""
InfiniteCylinder: Qsca vs Wavelength
====================================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of wavelength for cylindrical scatterers using PyMieSim, considering cylinders with different diameters and refractive indices.
"""
import numpy as np
from PyMieSim import (
    ureg,
    InfiniteCylinderSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)


polarization_state = PolarizationSet(
    angles=0 * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(400, 1000, 150)
    * ureg.nanometer,
    polarization=polarization_state,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=[200, 150] * ureg.nanometer,
    material=[2, 3, 4],
    medium=[1],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="source:wavelength", std="scatterer:material")
