"""
Sphere: A1 scattering coefficient
===================================

"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=400 * nanometer,
    polarization=[0] * degree,
    optical_power=1e-6 * watt,
    NA=0.2 * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(300, 1000, 100) * nanometer,
    property=[1.2, 1.25] * RIU,
    medium_property=[1.0] * RIU,
    source=source
)
#
# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
dataframe = experiment.get('a1')

# print(dataframe)
# %%
# Plotting the results
dataframe.plot_data(x='scatterer:diameter', std='scatterer:property')