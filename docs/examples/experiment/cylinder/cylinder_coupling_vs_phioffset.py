"""
InfiniteCylinder: Goniometer
============================

This example demonstrates how to use a goniometer setup to measure and visualize the coupling efficiency as a function of angular displacement for cylindrical scatterers using PyMieSim.
"""
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.detector_set import PhotodiodeSet
from PyMieSim.experiment.scatterer_set import InfiniteCylinderSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import print_available, SellmeierMaterial

print_available()

polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=[2000] * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.5, 0.3, 0.1, 0.05],
    phi_offset=np.linspace(-180, 180, 200) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[400],
    polarization_filter=None,
)

experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling", scale_unit=True)

dataframe.plot(x="detector:phi_offset")
