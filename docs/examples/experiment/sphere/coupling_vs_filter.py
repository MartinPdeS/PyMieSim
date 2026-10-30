"""
Sphere: Coupling vs polarization filter
=======================================
"""

import PyOptik
from PyMieSim import from_pyoptik
import numpy as np
from PyMieSim import (
    Measure,
    ureg,
    PhotodiodeSet,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)



polarization_set = PolarizationSet(
    angles=[0] * ureg.radian
)

source = GaussianSet(
    wavelength=[408] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=0.2,
)

scatterer = SphereSet(
    diameter=np.linspace(1000, 1200, 20) * ureg.nanometer,
    material=[from_pyoptik(PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7")), from_pyoptik(PyOptik.load_material(name="main/H2O/Daimon-19.0C"))],
    medium=[1],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.1],
    phi_offset=[-180] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    polarization_filter=np.linspace(-180, 180, 100) * ureg.degree,
    sampling=[500],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

result = experiment.run(Measure.COUPLING)

print("Result dimensions and shape:", result.dims, result.shape)
print("Measure units:", result.units)

result.plot(x="detector:polarization_filter", std="scatterer:diameter")
