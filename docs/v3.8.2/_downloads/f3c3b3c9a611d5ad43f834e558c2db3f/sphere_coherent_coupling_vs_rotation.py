"""
Sphere: Coherent mode field rotation
====================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from TypedUnit import ureg

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=1200 * ureg.nanometer,
    polarization=90 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=[2000, 2300] * ureg.nanometer,
    property=Material.BK7,
    medium_property=1 * ureg.RIU,
    source=source,
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number="HG11",
    NA=[0.05] * ureg.AU,
    phi_offset=0 * ureg.degree,
    gamma_offset=20 * ureg.degree,
    sampling=400 * ureg.AU,
    rotation=numpy.linspace(0, 180, 200) * ureg.degree,
    polarization_filter=None,
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get("coupling")

# %%
# Plotting the results
dataframe.plot(x="detector:rotation")
