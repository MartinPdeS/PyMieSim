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

polystyrene = SellmeierMaterial(material_name="BK7")

scatterer = SphereSet(
    diameter=np.linspace(10, 1000, 150) * ureg.nanometer,
    material=MaterialSet(materials=[polystyrene]),
    medium=MediumSet(refractive_indices=[1.33, 1.34, 1.5]),
)

experiment = Experiment(scatterer_set=scatterer, source_set=source)

result = experiment.get("Qsca")

result.plot(x="scatterer:diameter", show=True)
