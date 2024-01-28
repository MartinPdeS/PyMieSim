"""
Cylinder: Goniometer
====================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim.materials import BK7
from PyMieSim import measure

# %%
# Defining the source to be employed.
source_set = Gaussian(
    wavelength=1200e-9,
    polarization_value=90,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = Cylinder(
    diameter=2000e-9,
    material=BK7,
    n_medium=1,
    source_set=source_set
)

# %%
# Defining the detector to be employed.
detector_set = Photodiode(
    NA=[0.5, 0.3, 0.1, 0.05],
    phi_offset=np.linspace(-180, 180, 400),
    gamma_offset=0,
    sampling=400,
    polarization_filter=None
)

# %%
# Defining the experiment setup
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=detector_set
)

# %%
# Measuring the properties
data = experiment.get(measure.coupling)

# %%
# Plotting the results
figure = data.plot(
    x=detector_set.phi_offset,
    y_scale='log',
    normalize=True
)

_ = figure.show()
