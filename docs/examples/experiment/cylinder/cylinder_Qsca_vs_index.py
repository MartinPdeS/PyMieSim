"""
InfiniteCylinder: Qsca vs Index
===============================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of refractive index for cylindrical scatterers using PyMieSim, considering multiple wavelengths.
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
    wavelength=[500, 1000, 1500] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=[800] * ureg.nanometer,
    material=np.linspace(1.3, 1.9, 1500),
    medium=[1.0],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.get("Qsca", "Qext")

result.plot(x="scatterer:material")
