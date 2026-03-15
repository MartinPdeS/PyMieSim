"""
InfiniteCylinder: Qsca vs wavelength std
========================================

"""
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import InfiniteCylinderSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(200, 1800, 300) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)
scatterer = InfiniteCylinderSet(
    diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
    material=[1.4] * ureg.RIU,
    medium=[1.0] * ureg.RIU,
)

experiment = Setup(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="source:wavelength", std="scatterer:diameter")
