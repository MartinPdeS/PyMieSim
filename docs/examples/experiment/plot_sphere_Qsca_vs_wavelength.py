"""
Sphere: Qsca vs wavelength mean
===============================

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
    wavelength=np.linspace(400e-9, 1000e-9, 500),
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=[200e-9, 150e-9, 100e-9],
    index=[2, 4],
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
data = experiment.get(measure.Csca)

# %%
# Plotting the results
figure = data.plot(x=source.wavelength)

_ = figure.show()
