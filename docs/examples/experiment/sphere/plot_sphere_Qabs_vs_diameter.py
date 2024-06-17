"""
Sphere: Qsca vs diameter
========================

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
    wavelength=400e-9,
    polarization=0,
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(1e-09, 800e-9, 300),
    material=[UsualMaterial.Silver, UsualMaterial.Gold, UsualMaterial.Aluminium],
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
data = experiment.get(measure.Qabs)

# %%
# Plotting the results
figure = data.plot(
    x=scatterer.diameter,
    y_scale="log"
)

_ = figure.show()
