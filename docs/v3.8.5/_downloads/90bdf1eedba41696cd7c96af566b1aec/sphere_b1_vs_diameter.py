"""
Sphere: B1 scattering coefficient
=================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

source = Gaussian(
    wavelength=400 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)
scatterer = Sphere(
    diameter=np.linspace(100, 10000, 800) * ureg.nanometer,
    refractive_index=1.4 * ureg.RIU,
    medium_refractive_index=1 * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("b1")

dataframe.plot(x="scatterer:diameter")
