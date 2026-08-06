"""
Sphere: Goniometer
==================

"""
from PyMieSim import (
    ureg,
    PhotodiodeSet,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
    print_available,
    SellmeierMaterial,
)
import numpy


print_available()

polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[600, 1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)
scatterer = SphereSet(
    diameter=[1000] * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = PhotodiodeSet(
    numerical_aperture=[0.1, 0.15, 0.2],
    phi_offset=numpy.linspace(-180, 180, 200) * ureg.degree,
    cache_numerical_aperture=[0.05],
    gamma_offset=[0, 40] * ureg.degree,
    sampling=[400]
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

result = experiment.get("coupling")


result.plot(x="detector:phi_offset", y='coupling', std='detector:NA')
