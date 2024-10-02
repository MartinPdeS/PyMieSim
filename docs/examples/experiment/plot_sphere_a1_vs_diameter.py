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
from PyMieSim.experiment import measure
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
    diameter=np.linspace(100, 10000, 100) * nanometer,
    index=[1.2, 1.25] * RIU,
    medium_index=[0.9, 1.1] * RIU,
    source=source
)
#
# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
dataframe = experiment.get('Qsca', 'Qpr')
# print(dataframe)
# %%
# Plotting the results
dataframe.plot_data(x='diameter', std='index')