"""
Cylinder: Qsca vs wavelength
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure

# %%
# Defining the source to be employed.
# The source is always a plane wave in the LMT framework.
# Here we want to study differents wavelength.
# The amplitude is set to one per default.
source_set = Gaussian(
    wavelength=np.linspace(400e-9, 1000e-9, 500),
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = Cylinder(
    diameter=[200e-9, 150e-9, 100e-9],
    index=[2, 3, 4],
    n_medium=1,
    source_set=source_set
)


# %%
# Defining the experiment setup
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the properties
data = experiment.get(measure.Qsca)

data = data.mean(scatterer_set.index)

# %%
# Plotting the results
figure = data.plot(
    x=source_set.wavelength
)

_ = figure.show()
