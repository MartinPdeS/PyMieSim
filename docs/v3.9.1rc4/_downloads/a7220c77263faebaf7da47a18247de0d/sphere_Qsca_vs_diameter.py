"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian

from PyMieSim.experiment import Setup
from PyOptik import Material

Material.print_available()

source = Gaussian(
    wavelength=[405, 810] * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)

scatterer = Sphere(
    diameter=np.linspace(10, 1000, 150) * ureg.nanometer,
    medium_property=[1.33, 1.34, 1.5] * ureg.RIU,
    property=Material.polystyren,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca", scale_unit=True)

dataframe.plot(x="scatterer:diameter", show=True)
