"""
CoreShell: B1 vs Core Diameter
==============================

This example demonstrates how to compute and visualize the B1 scattering parameter as a function of core diameter for CoreShell scatterers using PyMieSim.
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
# The source is always a plane wave in the LMT framework.
# The amplitude is set to one per default.
source = Gaussian(
    wavelength=800 * nanometer,  # 800 nm
    polarization=0 * degree,  # Linear polarization angle in radians
    optical_power=1e-3 * watt,  # 1 milliwatt
    NA=0.2 * AU  # Numerical Aperture
)

# %%
# Defining the ranging parameters for the scatterer distribution
# Here we look at core/shell scatterers and use constant shell diameter
# with variable core diameter
scatterer = CoreShell(
    core_diameter=np.geomspace(100, 3000, 500) * nanometer,  # Geometrically spaced core diameters
    shell_width=800 * nanometer,  # Shell width of 800 nm
    core_property=1.6 * RIU,  # Refractive index of the core
    shell_property=Material.BK7,  # BK7 glass material for the shell
    medium_property=1 * RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Defining the experiment setup
# Integrating the defined source and scatterers into a single experimental setup.
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the B1 scattering parameter
dataframe = experiment.get('b1')

# %%
# Plotting the results
# Visualizing how the B1 parameter varies with the core diameter.
dataframe.plot_data(x="scatterer:core_diameter")