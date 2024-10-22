"""
Sphere: Coherent Goniometer
===========================

"""


# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=1200 * nanometer,
    polarization=90 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=2000 * nanometer,
    property=Material.BK7,
    medium_property=1 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = CoherentMode(
    mode_number='LP11',
    NA=[0.5, 0.3, 0.1, 0.05] * AU,
    phi_offset=numpy.linspace(-180, 180, 300) * degree,
    gamma_offset=0 * degree,
    sampling=400 * AU,
    polarization_filter=10 * degree,
    rotation=0 * degree,  # Rotation of the mode field
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the properties
dataframe = experiment.get('coupling')

# # %%
# # Plotting the results
dataframe.plot_data(x="detector:phi_offset")
