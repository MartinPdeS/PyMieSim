"""
Cylinder: Coupling vs wavelength
================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.detector import LPMode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup


from PyMieSim import measure
from PyMieSim.materials import BK7

# %%
# Defining the source to be employed.
source_set = Gaussian(
    wavelength=np.linspace(950e-9, 1050e-9, 300),
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the ranging parameters for the scatterer distribution
# Here we look at cylinder scatterers a set diameter, refractive index and medium.
scatterer_set = Cylinder(
    diameter=np.linspace(100e-9, 8000e-9, 5),
    material=BK7,
    n_medium=1,
    source_set=source_set
)

# %%
# Defining the detector to be employed.
detector_set = LPMode(
    mode_number="LP11",
    NA=[0.05, 0.01],
    phi_offset=-180,
    gamma_offset=0,
    polarization_filter=None,
    sampling=300
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
    x=source_set.wavelength,
    std=scatterer_set.diameter
)

_ = figure.show()
