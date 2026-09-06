"""
InfiniteCylinder: Qsca vs wavelength std
========================================

"""
import numpy as np
from PyMieSim import (
    ureg,
    InfiniteCylinderSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)


polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(200, 1800, 300) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1] * ureg.watt,
    numerical_aperture=[0.2],
)
scatterer = InfiniteCylinderSet(
    diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
    material=[1.4],
    medium=[1.0],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.get("Qsca")

result.plot(x="source:wavelength", std="scatterer:diameter")
