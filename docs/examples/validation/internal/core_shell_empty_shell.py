"""
Effect of Shell dimensions in equivalent medium
===============================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import degree, watt, AU, RIU, micrometer, nanometer

# Setup parameters
scatterer_diameter = 0.3 * micrometer  # Diameter of the scatterer in meters
scatterer_index = 1.4 * RIU  # Refractive index of the scatterer
source_wavelength = 1.2 * micrometer  # Wavelength of the source in meters

# Experiment source and scatterer setup
source = Gaussian(
    wavelength=1.2 * micrometer,
    polarization=[0, 90] * degree,
    optical_power=1 * watt,
    NA=0.2 * AU
)

scatterer = CoreShell(
    core_diameter=300 * nanometer,
    shell_thickness=np.linspace(100, 300, 100) * nanometer,
    core_property=1.4 * RIU,
    shell_property=1.3 * RIU,
    medium_property=1.3 * RIU,
    source=source
)


# Configure experiment
experiment = Setup(scatterer=scatterer, source=source)

# Gather data
# %%
dataframe = experiment.get('Csca')


ax = dataframe.plot_data(x='scatterer:shell_thickness', show=False)

ax.set_ylim([0, 4.0e-16])

plt.show()


# %%
# As it is supposed the scattering Cross-section should not be affected by
# the shell thickness as it's refractive index is same as the surrounding medium.
dataframe = experiment.get('Qsca')

ax = dataframe.plot_data(x='scatterer:shell_thickness', show=False)

plt.show()
# Similarly the scattering decrease as the scatterer becomes technically larger but the effects of the shell is for no account.
