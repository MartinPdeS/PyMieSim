"""
InfiniteCylinder: Qabs vs Diameter
==================================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim.
"""

import PyOptik
from PyMieSim import from_pyoptik
import numpy as np
from PyMieSim import (
    Measure,
    ureg,
    InfiniteCylinderSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)




polarization_set = PolarizationSet(
    angles=[30.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = InfiniteCylinderSet(
    diameter=np.linspace(1, 800, 300) * ureg.nanometer,
    material=[from_pyoptik(PyOptik.load_material(name="main/Ag/Johnson")), from_pyoptik(PyOptik.load_material(name="main/Au/Olmon-ev")), from_pyoptik(PyOptik.load_material(name="main/Al/Rakic"))],
    medium=[1],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.run(Measure.QABS, Measure.QSCA, Measure.QEXT)

print("Result dimensions and shape:", result["Qabs"].dims, result["Qabs"].shape)
print("Measure units:", {name: value.units for name, value in result.items()})

result["Qabs"].plot(x="scatterer:diameter", y="Qabs")
