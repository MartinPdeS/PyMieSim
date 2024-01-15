"""
Sphere: Goniometer
==================

"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment import SphereSet, SourceSet, Setup, PhotodiodeSet
from PyMieSim.materials import BK7
from PyMieSim import measure

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer_set = SphereSet(
    diameter=2000e-9,
    material=BK7,
    n_medium=1
)

# %%
# Defining the source to be employed.
source_set = SourceSet(
    wavelength=1200e-9,
    linear_polarization=90,
    optical_power=1e-3,
    NA=0.2
)

# %%
# Defining the detector to be employed.
detector_set = PhotodiodeSet(
    NA=[0.5, 0.3, 0.1, 0.05],
    phi_offset=numpy.linspace(-180, 180, 400),
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
data = experiment.Get(measure.coupling)

# %%
# Plotting the results
figure = data.plot(
    x=detector_set.phi_offset,
    y_scale='log',
    normalize=True
)

_ = figure.show()
