"""
Sphere: B1 scattering coefficient
=================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import SphereSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)
scatterer = SphereSet(
    diameter=np.linspace(100, 10000, 800) * ureg.nanometer,
    refractive_index=[1.4] * ureg.RIU,
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("b1")

dataframe.plot(x="scatterer:diameter")
