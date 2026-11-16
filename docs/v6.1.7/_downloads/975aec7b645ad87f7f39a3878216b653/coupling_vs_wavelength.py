"""
Sphere: Coupling vs wavelength
==============================
"""

import PyOptik
from PyMieSim import from_pyoptik
import numpy as np
from PyMieSim import (
    Measure,
    ureg,
    CoherentModeSet,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)



polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=np.linspace(950, 1050, 200) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = SphereSet(
    diameter=np.linspace(1400, 1500, 10) * ureg.nanometer,
    material=[from_pyoptik(PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7"))],
    medium=[1],
)

detector = CoherentModeSet(
    mode_number=["LP11"],
    numerical_aperture=[0.05, 0.03],
    phi_offset=[-180] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    polarization_filter=[0, 90] * ureg.degree,
    rotation=[0] * ureg.degree,
    sampling=[300],
    medium=[from_pyoptik(PyOptik.load_material(name="main/H2O/Daimon-19.0C"), medium=True)]
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

result = experiment.run(Measure.COUPLING)

print("Result dimensions and shape:", result.dims, result.shape)
print("Measure units:", result.units)

result.plot(x="source:wavelength", std="scatterer:diameter")
