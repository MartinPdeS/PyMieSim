"""
Effect of Shell dimensions in equivalent medium
===============================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
from PyMieSim.units import ureg

# PyMieSim imports
from PyMieSim.experiment.scatterer_set import CoreShellSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup

# Setup parameters
scatterer_diameter = 0.3 * ureg.micrometer  # Diameter of the scatterer in meters
scatterer_index = 1.4  # Refractive index of the scatterer
source_wavelength = 1.2 * ureg.micrometer  # Wavelength of the source in meters

polarization_set = PolarizationSet(
    angles=[0, 90] * ureg.degree,
)

# Experiment source and scatterer setup
source = GaussianSet(
    wavelength=[1200] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1] * ureg.watt,
    numerical_aperture=[0.2],
)

scatterer = CoreShellSet(
    core_diameter=[300] * ureg.nanometer,
    shell_thickness=np.linspace(100, 300, 100) * ureg.nanometer,
    core_material=[1.4],
    shell_material=[1.3],
    medium=[1.3],
)


experiment = Setup(scatterer_set=scatterer, source_set=source)

dataframe = experiment.get("Csca")

figure = dataframe.plot(x="scatterer:shell_thickness")


# %%
# As it is supposed the scattering Cross-section should not be affected by
# the shell thickness as it's refractive index is same as the surrounding medium.
dataframe = experiment.get("Qsca")

figure = dataframe.plot(x="scatterer:shell_thickness")

# Similarly the scattering decrease as the scatterer becomes technically larger but the effects of the shell is for no account.
