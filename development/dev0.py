"""
Sphere: Coupling vs diameter
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=488 * nanometer,
    polarization=90 * degree,
    optical_power=1e-3 * watt,
    NA=[0.1] * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=numpy.linspace(10, 200, 600) * nanometer,
    property=1.40 * RIU,
    medium_property=1.33 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = Photodiode(
    NA=1.2 * AU,
    phi_offset=[0.0, 90] * degree,
    gamma_offset=0.0 * degree,
    sampling=600 * AU,
    polarization_filter=None
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get('coupling', drop_unique_level=True, scale_unit=True)

# %%
# Plotting the results
dataframe.plot_data(x='scatterer:diameter')
