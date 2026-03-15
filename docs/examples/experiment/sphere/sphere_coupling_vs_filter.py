"""
Sphere: Coupling vs polarization filter
=======================================
"""
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.detector_set import PhotodiodeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import print_available, SellmeierMaterial

print_available()

polarization_set = PolarizationSet(
    angles=[0] * ureg.radian
)

source = GaussianSet(
    wavelength=[408] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=0.2 * ureg.AU,
)

scatterer = SphereSet(
    diameter=np.linspace(1000, 1200, 20) * ureg.nanometer,
    material=[SellmeierMaterial("BK7"), SellmeierMaterial("water")],
    medium=[1] * ureg.RIU,
)

detector = PhotodiodeSet(
    numerical_aperture=[0.1] * ureg.AU,
    phi_offset=[-180] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    polarization_filter=np.linspace(-180, 180, 100) * ureg.degree,
    sampling=[500],
)

experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:polarization_filter", std="scatterer:diameter")
