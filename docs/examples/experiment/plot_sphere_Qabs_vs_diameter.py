"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

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
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=np.linspace(1, 800, 300) * nanometer,
    material=[Material.silver, Material.gold, Material.aluminium],
    medium_index=1 * RIU,
    source=source
)

# %%
# Defining the experiment setup
experiment = Setup(
    scatterer=scatterer,
    source=source
)

# %%
# Measuring the properties
dataframe = experiment.get(measure.Qabs)

# %%
# Plotting the results
dataframe.plot_data(x='diameter')