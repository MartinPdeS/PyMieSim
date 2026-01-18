"""
Sphere: A1 scattering coefficient
===================================

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
    polarization=[0] * ureg.degree,
    optical_power=1e-6 * ureg.watt,
    NA=0.2 * ureg.AU,
)
scatterer = Sphere(
    diameter=np.linspace(300, 1000, 100) * ureg.nanometer,
    property=[1.2, 1.25] * ureg.RIU,
    medium_property=[1.0] * ureg.RIU,
    source=source,
)
#
experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("a1")

dataframe.plot(x="scatterer:diameter", std="scatterer:property")
