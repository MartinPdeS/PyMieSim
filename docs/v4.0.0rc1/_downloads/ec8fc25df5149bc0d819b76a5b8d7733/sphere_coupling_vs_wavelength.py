"""
Sphere: Coupling vs wavelength
==============================
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

source = Gaussian(
    wavelength=np.linspace(950, 1050, 200) * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)

scatterer = Sphere(
    diameter=np.linspace(100, 8000, 5) * ureg.nanometer,
    refractive_index=Material.BK7,
    medium_refractive_index=1 * ureg.RIU,
    source=source,
)

detector = CoherentMode(
    mode_number="LP11",
    NA=[0.05, 0.01] * ureg.AU,
    phi_offset=-180 * ureg.degree,
    gamma_offset=0 * ureg.degree,
    polarization_filter=[0, 90] * ureg.degree,
    rotation=0 * ureg.degree,
    sampling=300 * ureg.AU,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling", scale_unit=True)

dataframe.plot(x="source:wavelength", std="scatterer:diameter")
