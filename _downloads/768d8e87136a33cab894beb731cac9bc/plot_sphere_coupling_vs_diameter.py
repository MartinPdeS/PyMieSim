"""
Sphere: Coupling vs diameter
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim import measure
from PyMieSim.materials import BK7

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
scatterer_set = Sphere(
    diameter=numpy.linspace(100e-9, 3000e-9, 600),
    material=BK7,
    n_medium=1.0,
    source_set=source_set
)

# %%
# Defining the detector to be employed.
detector_set = Photodiode(
    NA=[0.15, 0.1, 0.05],
    phi_offset=-180.0,
    gamma_offset=0.0,
    sampling=600,
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
    x=scatterer_set.diameter,
    y_scale='linear',
    normalize=True
)

_ = figure.show()
