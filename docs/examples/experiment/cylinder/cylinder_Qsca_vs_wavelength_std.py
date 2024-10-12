"""
Cylinder: Qsca vs wavelength std
================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=np.linspace(200, 1800, 300) * nanometer,
    polarization=0 * degree,
    optical_power=1 * watt,
    NA=0.2 * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Cylinder(
    diameter=np.linspace(400, 1400, 10) * nanometer,
    property=Material.silver,
    medium_property=1 * RIU,
    source=source
)


# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
dataframe = experiment.get('Qsca')

# %%
# Plotting the results
dataframe.plot_data(x="source:wavelength", std='scatterer:diameter')

