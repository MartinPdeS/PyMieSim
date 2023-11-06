"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment import SphereSet, SourceSet, Setup
from PyMieSim import measure

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = SphereSet(
    diameter=np.geomspace(6.36e-09, 10000e-9, 200),
    index=1.4,
    n_medium=1
)

# %%
# Defining the source to be employed.
source_set = SourceSet(
    wavelength=[500e-9, 1000e-9, 1500e-9],
    linear_polarization=30,
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

# %%
# Plotting the results
figure = data.plot(
    x=scatterer_set.diameter
)

_ = figure.show()
