"""
Sphere: Qsca vs index
=====================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=[500e-9, 1000e-9, 1500e-9],
    polarization=30,
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=800e-9,
    index=np.linspace(1.3, 1.9, 1500),
    medium_index=1,
    source=source
)

# %%
# Defining the experiment setup
experiment = Setup(
    scatterer=scatterer,
    source=source
)

# %%
# Measuring the properties
data = experiment.get(measure.Qsca)

# %%
# Plotting the results
figure = data.plot(x=scatterer.index)

_ = figure.show()
