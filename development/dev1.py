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
from PyMieSim.polarization import RightCircular
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU, micrometer

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=550 * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    # diameter=[0.050, 0.473] * micrometer,
    diameter=numpy.linspace(50, 474, 5) * nanometer,
    property=(1.33 + 0.14) * RIU,
    medium_property=1.33 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = Photodiode(
    NA=[0.05] * AU,
    phi_offset=numpy.linspace(0, 180, 200) * degree,
    cache_NA=0.0 * AU,
    gamma_offset=0 * degree,
    sampling=400 * AU,
    polarization_filter=None
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get('coupling')

# %%
# Plotting the results
ax = dataframe.plot(x="detector:phi_offset", yscale='log')
