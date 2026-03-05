"""
Sphere: coherent coupling vs sampling
=====================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.units import ureg
from PyMieSim.experiment.detector import CoherentModeSet
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
    diameter=[100] * ureg.nanometer,
    material=[Material.BK7],
    medium_refractive_index=[1] * ureg.RIU,
    source=source,
)

detector = CoherentModeSet(
    mode_number=["LP01"],
    rotation=[0] * ureg.degree,
    numerical_aperture=[0.1] * ureg.AU,
    phi_offset=numpy.linspace(-20, 20, 200) * ureg.degree,
    gamma_offset=[0] * ureg.degree,
    sampling=[10, 20, 40, 80, 160, 500],
    polarization_filter=[0] * ureg.degree,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:phi_offset")
