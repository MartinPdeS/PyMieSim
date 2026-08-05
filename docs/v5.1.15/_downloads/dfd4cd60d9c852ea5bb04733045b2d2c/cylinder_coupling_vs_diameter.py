"""
InfiniteCylinder: Coupling vs Diameter
======================================

This example demonstrates how to compute and visualize the coupling efficiency as a function of diameter for cylindrical scatterers using PyMieSim.
"""

import numpy as np
from PyMieSim import (
    ureg,
    PhotodiodeSet,
    InfiniteCylinderSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)


polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[100, 1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(100, 300, 200) * ureg.nanometer,
    material=[1.4],
    medium=[1.0],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.1],
    phi_offset=[-180.0] * ureg.degree,
    gamma_offset=[0.0] * ureg.degree,
    sampling=[600],
    polarization_filter=None,
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="scatterer:diameter")
