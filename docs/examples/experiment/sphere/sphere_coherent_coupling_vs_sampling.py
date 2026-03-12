"""
Sphere: coherent coupling vs sampling
=====================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
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
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)
scatterer = SphereSet(
    diameter=[1000] * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1] * ureg.RIU,
)

detector = CoherentModeSet(
    mode_number=["LP01"],
    rotation=[0] * ureg.degree,
    numerical_aperture=[0.1] * ureg.AU,
    phi_offset=numpy.linspace(-80, 80, 200) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[10, 20, 40, 80, 160],
    polarization_filter=[0] * ureg.degree,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:phi_offset")
