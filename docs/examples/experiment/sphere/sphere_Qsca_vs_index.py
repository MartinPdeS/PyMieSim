"""
Sphere: Qsca vs index
=====================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=[500.0, 1000.0, 1500.0] * ureg.nanometer,
    polarization=30.0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=800.0 * ureg.nanometer,
    property=np.linspace(1.3, 1.9, 150) * ureg.RIU,
    medium_property=1.0 * ureg.RIU,
    source=source,
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
dataframe = experiment.get("Qsca")

# %%
# Plotting the results
dataframe.plot(x="scatterer:property")
