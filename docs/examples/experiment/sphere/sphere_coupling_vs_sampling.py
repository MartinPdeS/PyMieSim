"""
Sphere: Coupling vs sampling
============================

"""
import numpy
from PyMieSim import (
    ureg,
    PhotodiodeSet,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
    print_available,
    SellmeierMaterial,
    SellmeierMedium,
)



print_available()


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
    material=[SellmeierMaterial("BK7")],
    medium=[SellmeierMedium("water")],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.2],
    phi_offset=numpy.linspace(-20, 20, 400) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[20, 40, 80, 160],
    medium=[1.0],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

result = experiment.get("coupling")

result.plot(x="detector:phi_offset")
