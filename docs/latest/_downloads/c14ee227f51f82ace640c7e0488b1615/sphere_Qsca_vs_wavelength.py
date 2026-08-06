"""
Sphere: Qsca vs wavelength
==========================

"""
import numpy as np
from PyMieSim import (
    ureg,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
    print_available,
    SellmeierMaterial,
)


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

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.get("Qsca", "Qpr")

result.isel({"measure": 0}).plot(x="source:wavelength", y="Qsca")
