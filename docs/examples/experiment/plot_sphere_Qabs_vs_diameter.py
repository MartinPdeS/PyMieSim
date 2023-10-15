"""
Qsca vs diameter
================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment import SphereSet, SourceSet, Setup
from PyMieSim.materials import Gold, Silver, Aluminium
from PyMieSim import measure

# %%
# Defining the ranging parameters for the scatterer distribution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scatterer_set = SphereSet(
    diameter=np.linspace(1e-09, 800e-9, 300),
    material=[Silver, Gold, Aluminium],
    n_medium=1
)

# %%
# Defining the source to be employed.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source_set = SourceSet(
    wavelength=400e-9,
    polarization=0,
    amplitude=1
)

# %%
# Defining the experiment setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the properties
# ~~~~~~~~~~~~~~~~~~~~~~~~
data = experiment.Get(measures=measure.Qabs)

# %%
# Plotting the results
# ~~~~~~~~~~~~~~~~~~~~
figure = data.plot(
    y=measure.Qabs,
    x=scatterer_set.diameter,
    y_scale="log"
)

_ = figure.show()

# -
