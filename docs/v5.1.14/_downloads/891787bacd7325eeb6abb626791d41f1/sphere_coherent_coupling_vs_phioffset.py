"""
Sphere: Coherent Goniometer
===========================

"""
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.detector_set import CoherentModeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import print_available, SellmeierMaterial

print_available()

polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)
scatterer = SphereSet(
    diameter=[2000] * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = CoherentModeSet(
    mode_number="LP11",
    numerical_aperture=[0.5, 0.3, 0.1, 0.05],
    phi_offset=numpy.linspace(-180, 180, 300) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[400],
    polarization_filter=[10] * ureg.degree,
    rotation=[0] * ureg.degree,
)

experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:phi_offset")
