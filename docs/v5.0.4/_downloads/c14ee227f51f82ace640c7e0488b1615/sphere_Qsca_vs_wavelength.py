"""
Sphere: Qsca vs wavelength
==========================

"""
import numpy as np
from PyMieSim.units import ureg

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
    wavelength=np.linspace(400, 1000, 50) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

bk7 = SellmeierMaterial("BK7")

scatterer = SphereSet(
    diameter=[200] * ureg.nanometer,
    material=[bk7],
    medium=[1],
)

experiment = Setup(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca", "Qpr", scale_unit=True)

dataframe.plot(x="source:wavelength")
