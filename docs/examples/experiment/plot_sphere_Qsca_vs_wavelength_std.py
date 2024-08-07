"""
Sphere: Qsca vs wavelength STD
==============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import UsualMaterial
from PyMieSim.experiment import measure

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=np.linspace(200e-9, 1800e-9, 300),
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(400e-9, 1400e-9, 10),
    material=UsualMaterial.Silver,
    medium_index=1,
    source=source
)


# %%
# Defining the experiment setup
experiment = Setup(
    scatterer=scatterer,
    source=source
)

# %%
# Measuring the properties
data = experiment.get(measure.Qsca)

# %%
# Plotting the results
figure = data.plot(
    x=source.wavelength,
    y_scale='log',
    std=scatterer.diameter
)

_ = figure.show()

# -
