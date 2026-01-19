"""
Sphere: Qsca vs wavelength
==========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import MaterialBank

source = Gaussian(
    wavelength=np.linspace(400, 1000, 50) * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)

scatterer = Sphere(
    diameter=[200] * ureg.nanometer,
    refractive_index=MaterialBank.BK7,
    medium_refractive_index=1 * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca", "Qpr", scale_unit=True)

dataframe.plot(x="source:wavelength")
