"""
Sphere: Coupling vs sampling
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
from PyMieSim.experiment import measure
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=400 * nanometer,
    polarization=90 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=5000 * nanometer,
    material=Material.BK7,
    medium_index=1 * RIU,
    source=source
)

# %%
# Defining the detector to be employed.
detector = Photodiode(
    NA=[0.2] * AU,
    phi_offset=numpy.linspace(-20, 20, 400) * degree,
    gamma_offset=0 * degree,
    sampling=[20, 40, 80, 160] * AU,
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
dataframe = experiment.get(measure.coupling)

# %%
# Plotting the results
dataframe.plot_data(x="phi_offset")
