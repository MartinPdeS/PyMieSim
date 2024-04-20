"""
Sphere: Goniometer
==================

"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import UsualMaterial
from PyMieSim import measure

# %%
# Defining the source to be employed.
source_set = Gaussian(
    wavelength=635e-9,
    polarization_value=[0, 90],
    polarization_type='linear',
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = Sphere(
    diameter=numpy.linspace(1000e-9 - 100e-9, 1000e-9 + 100e-9, 10),
    material=UsualMaterial.Polystyrene,
    medium_material=UsualMaterial.Water,
    source_set=source_set
)

# %%
# Defining the detector to be employed.
detector_set = Photodiode(
    NA=[0.1, 0.01],
    phi_offset=numpy.linspace(-165, 165, 200),
    gamma_offset=0,
    sampling=1000,
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
    x=experiment.phi_offset,
    y_scale='log',
    std=experiment.diameter
    # normalize=True
)

_ = figure.show()
