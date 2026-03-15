"""
Sphere: Qsca vs diameter
========================

"""
import numpy as np

from PyMieSim.units import ureg
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.material_set import MaterialSet, MediumSet
from PyMieSim.experiment.setup import Setup
from PyMieSim.material import print_available, SellmeierMaterial

print_available()

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[405] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

polystyrene = SellmeierMaterial("polystyrene")

scatterer = SphereSet(
    diameter=np.linspace(10, 1000, 150) * ureg.nanometer,
    material=MaterialSet([polystyrene]),
    medium=MediumSet([1.33, 1.34, 1.5] * ureg.RIU),
)

experiment = Setup(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca", scale_unit=True)

dataframe.plot(x="scatterer:diameter", show=True)
