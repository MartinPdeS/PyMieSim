"""
Sphere: Coupling vs polarization filter
=======================================
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
    wavelength=950e-9,
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(100e-9, 2000e-9, 20),
    material=UsualMaterial.BK7,
    medium_index=1,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number="HG11:00",
    NA=[0.1],
    phi_offset=-180,
    gamma_offset=0,
    polarization_filter=np.linspace(-180, 180, 100),
    sampling=300,
    rotation=0,  # Rotation of the mode field
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
    x=detector.polarization_filter,
    std=scatterer.diameter
)

_ = figure.show()
