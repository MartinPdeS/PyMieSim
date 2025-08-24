"""
Sphere: B1 scattering coefficient
=================================

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
    wavelength=400 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(100, 10000, 800) * ureg.nanometer,
    property=1.4 * ureg.RIU,
    medium_property=1 * ureg.RIU,
    source=source
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
dataframe = experiment.get('b1')

# %%
# Plotting the results
dataframe.plot(x='scatterer:diameter')
