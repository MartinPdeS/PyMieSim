"""
Sphere: Qsca vs wavelength
==========================

"""

import PyOptik
from PyMieSim import from_pyoptik
import numpy as np
from PyMieSim import (
    Measure,
    ureg,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)



polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(400, 1000, 50) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

bk7 = from_pyoptik(PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7"))

scatterer = SphereSet(
    diameter=[200] * ureg.nanometer,
    material=[bk7],
    medium=[1],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.run(Measure.QSCA, Measure.QPR)

print("Result dimensions and shape:", result["Qsca"].dims, result["Qsca"].shape)
print("Measure units:", {name: value.units for name, value in result.items()})

result["Qsca"].plot(x="source:wavelength", y="Qsca")
