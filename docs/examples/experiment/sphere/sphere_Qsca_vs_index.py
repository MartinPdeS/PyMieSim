"""
Sphere: Qsca vs index
=====================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

source = Gaussian(
    wavelength=[500.0, 1000.0, 1500.0] * ureg.nanometer,
    polarization=30.0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)

scatterer = Sphere(
    diameter=800.0 * ureg.nanometer,
    property=np.linspace(1.3, 1.9, 150) * ureg.RIU,
    medium_property=1.0 * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="scatterer:property")
