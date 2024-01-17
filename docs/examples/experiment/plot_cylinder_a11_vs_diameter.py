"""
Cylinder: A1 scattering coefficient
===================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = Cylinder(
    diameter=np.linspace(100e-9, 10000e-9, 800),
    index=1.4,
    n_medium=1
)

# %%
# Defining the source to be employed.
source_set = Gaussian(
    wavelength=400e-9,
    linear_polarization=90,
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the experiment setup
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the properties
data = experiment.Get(measure.a21)

# %%
# Plotting the results
figure = data.plot(
    x=scatterer_set.diameter
)

_ = figure.show()
