"""
Sphere: Coupling vs sampling
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.detector_set import PhotodiodeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import print_available, SellmeierMaterial, SellmeierMedium

print_available()


polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)
scatterer = SphereSet(
    diameter=[5000] * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[SellmeierMedium("water")],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.2] * ureg.AU,
    phi_offset=numpy.linspace(-20, 20, 400) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[20, 40, 80, 160]
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:phi_offset")
