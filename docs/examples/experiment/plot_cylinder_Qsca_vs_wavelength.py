"""
Cylinder: Qsca vs wavelength
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment import CylinderSet, SourceSet, Setup
from PyMieSim import measure

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = CylinderSet(
    diameter=[200e-9, 150e-9, 100e-9],
    index=[2, 3, 4],
    n_medium=1
)

# %%
# Defining the source to be employed.
# The source is always a plane wave in the LMT framework.
# Here we want to study differents wavelength.
# The amplitude is set to one per default.
source_set = SourceSet(
    wavelength=np.linspace(400e-9, 1000e-9, 500),
    linear_polarization=0,
    amplitude=1
)

# %%
# Defining the experiment setup
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the properties
data = experiment.Get(measure.Qsca)

data = data.mean(scatterer_set.index)

# %%
# Plotting the results
figure = data.plot(
    x=source_set.wavelength
)

_ = figure.show()
