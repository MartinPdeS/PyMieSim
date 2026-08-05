"""
Sphere: coherent coupling vs sampling
=====================================

"""
import numpy
from PyMieSim import (
    ureg,
    CoherentModeSet,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
    print_available,
    SellmeierMaterial,
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
    diameter=[1000] * ureg.nanometer,
    material=[SellmeierMaterial("BK7")],
    medium=[1],
)

detector = CoherentModeSet(
    mode_number=["LP01"],
    rotation=[0] * ureg.degree,
    numerical_aperture=[0.1],
    phi_offset=numpy.linspace(-80, 80, 200) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[10, 20, 40, 80, 160],
    polarization_filter=[0] * ureg.degree,
)

experiment = Experiment(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:phi_offset")
