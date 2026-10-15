"""
Sphere: Coupling vs sampling
============================

"""

import PyOptik
from PyMieSim import from_pyoptik
import numpy
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
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = SphereSet(
    diameter=[5000] * ureg.nanometer,
    material=[from_pyoptik(PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7"))],
    medium=[from_pyoptik(PyOptik.load_material(name="main/H2O/Daimon-19.0C"), medium=True)],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.2],
    phi_offset=numpy.linspace(-20, 20, 400) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[20, 40, 80, 160],
    medium=[1.0],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

result = experiment.get(Measure.COUPLING)

print("Result dimensions and shape:", result.dims, result.shape)
print("Measure units:", result.attrs["units"])

result.plot(x="detector:phi_offset")
