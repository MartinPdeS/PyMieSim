"""
Sphere: Qsca vs diameter
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
    MaterialSet,
    MediumSet,
    Experiment,
)


polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[405] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

bk7 = from_pyoptik(PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7"))

scatterer = SphereSet(
    diameter=np.linspace(10, 1000, 150) * ureg.nanometer,
    material=MaterialSet(materials=[bk7]),
    medium=MediumSet(refractive_indices=[1.33, 1.34, 1.5]),
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.get(Measure.QSCA)

print("Result dimensions and shape:", result.dims, result.shape)
print("Measure units:", result.attrs["units"])

result.plot(x="scatterer:diameter", show=True)
