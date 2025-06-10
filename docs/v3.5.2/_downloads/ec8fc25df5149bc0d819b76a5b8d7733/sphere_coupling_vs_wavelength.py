"""
Sphere: Coupling vs wavelength
==============================
"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=np.linspace(950, 1050, 200) * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)

# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(100, 8000, 5) * nanometer,
    property=Material.BK7,
    medium_property=1 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number="LP11",
    NA=[0.05, 0.01] * AU,
    phi_offset=-180 * degree,
    gamma_offset=0 * degree,
    polarization_filter=[0, 90] * degree,
    rotation=0 * degree,
    sampling=300 * AU
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get('coupling', scale_unit=True)

# %%
# Plotting the results
dataframe.plot(x="source:wavelength", std='scatterer:diameter')
