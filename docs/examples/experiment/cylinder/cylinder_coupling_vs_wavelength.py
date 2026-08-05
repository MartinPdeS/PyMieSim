"""
InfiniteCylinder: Coupling vs Wavelength
========================================

This example demonstrates how to compute and visualize the coupling efficiency as a function of wavelength for cylindrical scatterers using PyMieSim.
"""
import numpy as np
from PyMieSim import (
    ureg,
    CoherentModeSet,
    InfiniteCylinderSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
    print_available,
    SellmeierMaterial,
)


print_available()

polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(950, 1050, 300) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(300, 500, 5) * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = CoherentModeSet(
    mode_number=["LP11"],
    numerical_aperture=[0.05, 0.01],
    phi_offset=[-180] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[300],
    rotation=[0] * ureg.degree,
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling", scale_unit=True)

dataframe.plot(x="source:wavelength", std="scatterer:diameter")
