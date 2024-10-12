"""
Sphere: Qsca vs diameter
========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np

from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=[500] * nanometer,
    polarization=30 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = CoreShell(
    core_diameter=np.geomspace(6.36, 10000, 1500) * nanometer,
    shell_width=200 * nanometer,
    medium_property=1 * RIU,
    core_property=1.6 * RIU,
    shell_property=1.6 * RIU,
    source=source
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
dataframe = experiment.get('Csca', 'Cabs', scale_unit=True, drop_unique_level=True)

# %%
# Plotting the results
dataframe.plot_data(x='scatterer:core_diameter')