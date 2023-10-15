"""
Sphere: B1 scattering coefficient
=================================

"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment import SphereSet, SourceSet, Setup
from PyMieSim import measure

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = SphereSet(
    diameter=np.linspace(100e-9, 10000e-9, 800),
    index=1.4,
    n_medium=1
)

# %%
# Defining the source to be employed.
source_set = SourceSet(
    wavelength=400e-9,
    polarization=0,
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
data = experiment.Get(measures=[measure.b1])

# %%
# Plotting the results
figure = data.plot(
    y=measure.b1,
    x=scatterer_set.diameter
)

_ = figure.show()
