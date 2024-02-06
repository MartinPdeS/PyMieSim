"""
CoreShell: B1 vs core diameter
==============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.materials import BK7
from PyMieSim import measure

# %%
# Defining the source to be employed.
# The source is always a plane wave in the LMT framework.
# The amplitude is set to one per default.
source_set = Gaussian(
    wavelength=800e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the ranging parameters for the scatterer distribution
# Here we look at core/shell scatterers and use constant shell diameter
# with variable core diameter
scatterer_set = CoreShell(
    core_diameter=np.geomspace(100e-09, 3000e-9, 5000),
    shell_width=800e-9,
    core_index=1.6,
    shell_material=BK7,
    n_medium=1,
    source_set=source_set
)

# %%
# Defining the experiment setup
# With the source and scatterers defined we set them together
# in an experiment.
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the properties
# We are interesting here in the b_1 (first magnetic coefficient) parameter.
data = experiment.get(measure.a3)

# %%
# Plotting the results
figure = data.plot(
    x=scatterer_set.core_diameter,
    y_scale='linear'
)

_ = figure.show()
