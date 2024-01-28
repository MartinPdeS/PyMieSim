"""
Cylinder: Qsca vs diameter
==========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim.materials import Gold, Silver, Aluminium
from PyMieSim import measure

# %%
# Defining the source to be employed.
source_set = Gaussian(
    wavelength=400e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = Cylinder(
    diameter=np.linspace(1e-09, 800e-9, 300),
    material=[Silver, Gold, Aluminium],
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
data = experiment.get(measure.Qabs)

# %%
# Plotting the results
figure = data.plot(
    x=scatterer_set.diameter,
    y_scale="linear"
)

_ = figure.show()
