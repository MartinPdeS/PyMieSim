"""
Sphere: Qsca vs wavelength
==========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyOptik import MaterialBank
# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=np.linspace(400, 1000, 50) * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=[200] * nanometer,
    property=MaterialBank.BK7,
    medium_property=1 * RIU,
    source=source
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
dataframe = experiment.get('Qsca', scale_unit=True)

# %%
# Plotting the results
dataframe.plot_data(x="source:wavelength")
