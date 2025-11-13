"""
Sphere: Coupling vs polarization filter
=======================================
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=[950, 1050] * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=[1e-3] * ureg.watt,
    NA=0.2 * ureg.AU,
)


# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(100, 2000, 20) * ureg.nanometer,
    property=[Material.BK7, Material.water],
    medium_property=1 * ureg.RIU,
    source=source,
)

# %%
# Defining the detector to be employed.
detector = Photodiode(
    NA=[0.1] * ureg.AU,
    phi_offset=-180 * ureg.degree,
    gamma_offset=0 * ureg.degree,
    polarization_filter=np.linspace(-180, 180, 100) * ureg.degree,
    sampling=[500] * ureg.AU,
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get("coupling")

dataframe.plot(x="detector:polarization_filter", std="scatterer:diameter")
