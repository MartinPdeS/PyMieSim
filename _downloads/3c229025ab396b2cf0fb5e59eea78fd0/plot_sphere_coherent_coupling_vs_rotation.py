"""
Sphere: Coherent mode field rotation
====================================

"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import UsualMaterial
from PyMieSim.experiment import measure

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=1200e-9,
    polarization=90,
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=[2000e-9, 2300e-9],
    material=UsualMaterial.BK7,
    medium_index=1,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number='HG11',
    NA=[0.05],
    phi_offset=0,
    gamma_offset=20,
    sampling=400,
    rotation=numpy.linspace(0, 180, 200),
    polarization_filter=None
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
figure = data.plot(x=detector.rotation)

_ = figure.show()
