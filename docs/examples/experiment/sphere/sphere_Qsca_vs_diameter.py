"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian, PolarizationSet

from PyMieSim.experiment import Setup
from PyOptik import Material

import PyMieSim
PyMieSim.debug_mode = True

Material.print_available()

polarization_set = PolarizationSet(
    angles=[0.0] * ureg.degree,
)

source = Gaussian(
    wavelength=[405] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=1e-3 * ureg.watt,
    numerical_aperture=0.2 * ureg.AU,
)

scatterer = Sphere(
    diameter=np.linspace(10, 1000, 150) * ureg.nanometer,
    medium_refractive_index=[1.33, 1.34, 1.5] * ureg.RIU,
    refractive_index=Material.polystyren,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get("Qsca", scale_unit=True)

dataframe.plot(x="scatterer:diameter", show=True)
