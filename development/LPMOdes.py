"""
Sphere: Coupling vs wavelength
==============================
"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

from PyMieSim import measure
from PyMieSim.materials import BK7

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=632.8e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=1000e-9,
    index=1.4,
    medium_index=1.21,
    source=source
)

print(scatterer.print_properties())
