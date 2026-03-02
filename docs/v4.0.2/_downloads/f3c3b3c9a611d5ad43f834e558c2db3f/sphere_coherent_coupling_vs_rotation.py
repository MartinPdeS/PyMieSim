"""
Sphere: Coherent mode field rotation
====================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import Material

polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = Gaussian(
    wavelength=1200 * ureg.nanometer,
    polarization=polarization_set,
    optical_power=1e-3 * ureg.watt,
    numerical_aperture=0.2 * ureg.AU,
)
scatterer = Sphere(
    diameter=[2000, 2300] * ureg.nanometer,
    refractive_index=Material.BK7,
    medium_refractive_index=1 * ureg.RIU,
    source=source,
)

detector = CoherentMode(
    mode_number="HG11",
    numerical_aperture=[0.05] * ureg.AU,
    phi_offset=0 * ureg.degree,
    gamma_offset=20 * ureg.degree,
    sampling=400,
    rotation=numpy.linspace(0, 180, 200) * ureg.degree,
    polarization_filter=None,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:rotation")
