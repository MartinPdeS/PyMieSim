"""
Sphere: Coupling vs sampling
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import PhotodiodeSet
from PyMieSim.experiment.scatterer import SphereSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import Material


polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[400] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)
scatterer = SphereSet(
    diameter=[5000] * ureg.nanometer,
    material=[Material.BK7],
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

detector = PhotodiodeSet(
    numerical_aperture=[0.2] * ureg.AU,
    phi_offset=numpy.linspace(-20, 20, 400) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[20, 40, 80, 160]
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:phi_offset")
