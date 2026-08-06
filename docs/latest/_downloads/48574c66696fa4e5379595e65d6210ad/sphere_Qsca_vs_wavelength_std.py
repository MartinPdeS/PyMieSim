"""
Sphere: Qsca vs wavelength STD
==============================

"""
import numpy as np
from PyMieSim import (
    ureg,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
    print_available,
    TabulatedMaterial,
)


print_available()

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(200, 1800, 200) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

silver = TabulatedMaterial("silver")
gold = TabulatedMaterial("gold")

scatterer = SphereSet(
    diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
    material=[silver, gold, 1.4],
    medium=[1],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.get("Qsca")

result.plot(x="source:wavelength", std="scatterer:diameter")

