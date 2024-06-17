"""
Sphere: B1 scattering coefficient
=================================

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
    wavelength=400e-9,
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(100e-9, 10000e-9, 800),
    index=1.4,
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
data = experiment.get(measure.b1)

# %%
# Plotting the results
figure = data.plot(x=scatterer.diameter)

_ = figure.show()
