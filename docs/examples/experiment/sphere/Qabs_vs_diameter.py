"""
Sphere: Qabs vs diameter
========================

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
    angles=[0] * ureg.degree
)

source = GaussianSet(
    wavelength=[400, 700] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

silver = from_pyoptik(PyOptik.load_material(name="main/Ag/Johnson"))

scatterer = SphereSet(
    diameter=np.linspace(1, 800, 300) * ureg.nanometer,
    material=[silver],
    medium=[1],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.run(Measure.QSCA)

print("Result dimensions and shape:", result.dims, result.shape)
print("Measure units:", result.units)

result.plot(x="scatterer:diameter")
