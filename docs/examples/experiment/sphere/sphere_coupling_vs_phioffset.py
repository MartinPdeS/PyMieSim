"""
Sphere: Goniometer
==================

"""
from PyMieSim.units import ureg
import numpy

from PyMieSim.experiment.detector_set import PhotodiodeSet
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
    diameter=[20] * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.1, 0.2],
    phi_offset=numpy.linspace(-180, 180, 200) * ureg.degree,
    cache_numerical_aperture=[0.05],
    gamma_offset=[0] * ureg.degree,
    sampling=[400]
)

experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:phi_offset")
