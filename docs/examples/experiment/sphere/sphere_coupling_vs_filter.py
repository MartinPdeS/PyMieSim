"""
Sphere: Coupling vs polarization filter
=======================================
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import PhotodiodeSet
from PyMieSim.experiment.scatterer import SphereSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import Material

polarization_set = PolarizationSet(
    angles=[0] * ureg.radian
)

source = GaussianSet(
    wavelength=[950, 1050] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=0.2 * ureg.AU,
)

scatterer = SphereSet(
    diameter=np.linspace(100, 2000, 20) * ureg.nanometer,
    material=[Material.BK7, Material.water],
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

detector = PhotodiodeSet(
    numerical_aperture=[0.1] * ureg.AU,
    phi_offset=[-180] * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    polarization_filter=np.linspace(-180, 180, 100) * ureg.degree,
    sampling=[500],
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:polarization_filter", std="scatterer:diameter")
