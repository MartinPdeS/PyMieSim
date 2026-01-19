"""
Cylinder: Qsca vs wavelength std
================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

source = Gaussian(
    wavelength=np.linspace(200, 1800, 300) * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.2 * ureg.AU,
)
scatterer = Cylinder(
    diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
    refractive_index=Material.silver,
    medium_refractive_index=1 * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca")

dataframe.plot(x="source:wavelength", std="scatterer:diameter")
