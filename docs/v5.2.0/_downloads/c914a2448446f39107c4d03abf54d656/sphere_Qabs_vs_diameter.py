"""
Sphere: Qabs vs diameter
========================

"""
import numpy as np
from PyMieSim import (
    ureg,
    SphereSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
    print_available,
    TabulatedMaterial,
)


print_available()

polarization_set = PolarizationSet(
    angles=[0] * ureg.degree
)

source = GaussianSet(
    wavelength=[400, 700] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

silver = TabulatedMaterial("silver")

scatterer = SphereSet(
    diameter=np.linspace(1, 800, 300) * ureg.nanometer,
    material=[silver],
    medium=[1],
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="scatterer:diameter", yscale="log")
