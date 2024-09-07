"""
Sphere: Coupling vs sampling
============================

"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.experiment import measure

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=400e-9,
    polarization=90,
    optical_power=1e-3,
    NA=0.2
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=5000e-9,
    material=Material.BK7,
    medium_index=1,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number='LP01',
    rotation=0,
    NA=[0.2],
    phi_offset=numpy.linspace(-20, 20, 400),
    gamma_offset=0,
    sampling=[20, 40, 80, 160, 1000, 2000],
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
data.plot(x=detector.phi_offset)
