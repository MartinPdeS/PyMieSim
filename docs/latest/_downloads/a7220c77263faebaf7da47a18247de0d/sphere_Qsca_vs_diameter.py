"""
Sphere: Qsca vs diameter
========================

"""
import numpy as np

from PyMieSim import (
    ureg,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    MaterialSet,
    MediumSet,
    Experiment,
    print_available,
    SellmeierMaterial,
)

print_available()

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[405] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

polystyrene = SellmeierMaterial("BK7")

scatterer = SphereSet(
    diameter=np.linspace(10, 1000, 150) * ureg.nanometer,
    material=MaterialSet([polystyrene]),
    medium=MediumSet([1.33, 1.34, 1.5]),
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca", scale_unit=True)

dataframe.plot(x="scatterer:diameter", show=True)
