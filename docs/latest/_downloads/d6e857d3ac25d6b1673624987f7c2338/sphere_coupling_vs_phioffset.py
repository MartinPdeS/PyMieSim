"""
Sphere: Goniometer
==================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
from PyMieSim.units import ureg
import numpy

from PyMieSim.experiment.detector import Photodiode
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
    diameter=20 * ureg.nanometer,
    refractive_index=Material.BK7,
    medium_refractive_index=1 * ureg.RIU,
    source=source,
)

detector = Photodiode(
    numerical_aperture=[0.1, 0.2] * ureg.AU,
    phi_offset=numpy.linspace(-180, 180, 200) * ureg.degree,
    cache_numerical_aperture=0.05 * ureg.AU,
    gamma_offset=0 * ureg.degree,
    sampling=400
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling")

dataframe.plot(x="detector:phi_offset")
