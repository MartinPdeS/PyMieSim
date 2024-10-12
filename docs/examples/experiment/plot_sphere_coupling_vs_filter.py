"""
Sphere: Coupling vs polarization filter
=======================================
"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import degree, nanometer, RIU, watt, AU
from PyOptik import Material

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=[950, 1050] * nanometer,
    polarization=0 * degree,
    optical_power=[1e-3] * watt,
    NA=0.2 * AU
)


# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(100, 2000, 20) * nanometer,
    property=[Material.BK7, Material.water],
    medium_property=1 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = Photodiode(
    NA=[0.1] * AU,
    phi_offset=-180 * degree,
    gamma_offset=0 * degree,
    polarization_filter=np.linspace(-180, 180, 100) * degree,
    sampling=[500] * AU,
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get('coupling')

dataframe.plot_data(x='detector:polarization_filter', std='scatterer:diameter')
