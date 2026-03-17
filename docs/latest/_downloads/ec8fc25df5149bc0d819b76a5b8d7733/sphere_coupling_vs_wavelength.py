"""
Sphere: Coupling vs wavelength
==============================
"""
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.detector_set import CoherentModeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import print_available, SellmeierMaterial, SellmeierMedium

print_available()

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(950, 1050, 200) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = SphereSet(
    diameter=np.linspace(1400, 1500, 10) * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = CoherentModeSet(
    mode_number=["LP11"],
    numerical_aperture=[0.05, 0.03],
    phi_offset=[-180] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    polarization_filter=[0, 90] * ureg.degree,
    rotation=[0] * ureg.degree,
    sampling=[300],
    medium=[SellmeierMedium("water")]
)

experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling", scale_unit=True)

dataframe.plot(x="source:wavelength", std="scatterer:diameter")
