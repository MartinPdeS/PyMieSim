"""
Sphere: Goniometer
==================

"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import UsualMaterial
from PyMieSim import measure

# %%
# Defining the source to be employed.
source_set = Gaussian(
    wavelength=635e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = Sphere(
    diameter=numpy.linspace(100e-9, 10e-6, 200),
    index=numpy.linspace(1.4, 1.6, 4),
    medium_material=UsualMaterial.Water,
    source_set=source_set
)

# %%
# Defining the experiment setup
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=None
)

# %%
# Measuring the properties
data = experiment.get(measure.Qback)

# %%
# Plotting the results
figure = data.plot(
    x=experiment.diameter,
    # y_scale='log',
    # normalize=True
)

_ = figure.show()
