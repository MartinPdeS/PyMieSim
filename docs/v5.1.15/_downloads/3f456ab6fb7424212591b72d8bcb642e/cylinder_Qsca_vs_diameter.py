"""
InfiniteCylinder: Qsca vs Diameter
==================================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim, considering multiple wavelengths.
"""
import numpy as np
from PyMieSim import (
    ureg,
    InfiniteCylinderSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)


polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[500, 1000] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=np.geomspace(6.36, 10000, 1000) * ureg.nanometer,
    material=[1.4,],
    medium=[1.0],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="scatterer:diameter")
