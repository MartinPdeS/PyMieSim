"""
Sphere: Coupling vs diameter
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure
from PyOptik import UsualMaterial

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=1200e-9,
    polarization=90,
    optical_power=1e-3,
    NA=[0.1]
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=numpy.linspace(100e-9, 3000e-9, 600),
    material=UsualMaterial.BK7,
    medium_index=1.0,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number='LP11',
    NA=numpy.linspace(0.05, 0.2, 5),
    rotation=0,
    phi_offset=[-180.0],
    gamma_offset=[0.0],
    sampling=[600],
    polarization_filter=[None],
    mean_coupling=False
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
    x=scatterer.diameter,
    y_scale='linear',
)

_ = figure.show()
