"""
Sphere: Coherent mode field rotation
====================================

"""

import PyOptik
from PyMieSim import from_pyoptik
import numpy
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
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = SphereSet(
    diameter=[2000, 2300] * ureg.nanometer,
    material=[from_pyoptik(PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7"))],
    medium=[1],
)

detector = CoherentModeSet(
    mode_number="HG11",
    numerical_aperture=[0.05],
    phi_offset=[0] * ureg.degree,
    gamma_offset=[20] * ureg.degree,
    sampling=[400],
    rotation=numpy.linspace(0, 180, 200) * ureg.degree,
    polarization_filter=None,
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

result = experiment.get(Measure.COUPLING)

print("Result dimensions and shape:", result.dims, result.shape)
print("Measure units:", result.attrs["units"])

result.plot(x="detector:rotation")
