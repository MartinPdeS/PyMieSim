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
from PyMieSim.material import print_available, SellmeierMaterial

print_available()

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(950, 1050, 200) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

scatterer = SphereSet(
    diameter=np.linspace(100, 8000, 10) * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1] * ureg.RIU,
)

detector = CoherentModeSet(
    mode_number=["LP11"],
    numerical_aperture=[0.05, 0.01] * ureg.AU,
    phi_offset=[-180] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    polarization_filter=[0, 90] * ureg.degree,
    rotation=[0] * ureg.degree,
    sampling=[300],
    medium=[1] * ureg.RIU,
)

experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling", scale_unit=True)

dataframe.plot(x="source:wavelength", std="scatterer:diameter")
