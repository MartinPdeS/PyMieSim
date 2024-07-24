"""
Sphere: Coupling vs wavelength
==============================
"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure
from PyOptik import UsualMaterial

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=np.linspace(950e-9, 1050e-9, 200),
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(100e-9, 8000e-9, 5),
    material=UsualMaterial.BK7,
    medium_index=1,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number="LP11",
    NA=[0.05, 0.01],
    phi_offset=-180,
    gamma_offset=0,
    polarization_filter=[0, None],
    rotation=0,
    sampling=300
)

# %%
# Defining the experiment setup
experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

# %%
# Measuring the properties
data = experiment.get(measure.coupling)

# %%
# Plotting the results
figure = data.plot(
    x=source.wavelength,
    std=scatterer.diameter
)

_ = figure.show()
