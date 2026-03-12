"""
Sphere: Qabs vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import print_available, TabulatedMaterial

print_available()

polarization_set = PolarizationSet(
    angles=[0] * ureg.degree
)

source = GaussianSet(
    wavelength=[400, 700] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

silver = TabulatedMaterial("silver")

scatterer = SphereSet(
    diameter=np.linspace(1, 800, 300) * ureg.nanometer,
    material=[silver],
    medium=[1] * ureg.RIU,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="scatterer:diameter", yscale="log")
