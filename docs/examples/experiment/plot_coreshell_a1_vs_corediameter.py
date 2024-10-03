"""
CoreShell: A1 vs Core Diameter
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
# Defining the source
# In the LMT framework, the source is always considered a plane wave with a default amplitude of one.
source = Gaussian(
    wavelength=800 * nanometer,  # 800 nm
    polarization=0 * degree,  # Linear polarization angle in radians
    optical_power=1e-3 * watt,  # 1 milliwatt
    NA=0.2 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
# Here, we explore core/shell scatterers with a constant shell diameter and variable core diameter.
scatterer = CoreShell(
    core_diameter=np.geomspace(100, 300, 10) * nanometer,  # Geometrically spaced core diameters
    shell_width=800 * nanometer,  # Shell width of 800 nm
    core_property=[1.4] * RIU,  # Refractive index of the core
    shell_property=[Material.BK7, Material.silver],  # BK7 glass material for the shell
    medium_property=1 * RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Defining the experiment setup
# Integrating the defined source and scatterers into a single experimental setup.
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the B1 scattering parameter
# Here, we're interested in the a3 (first magnetic coefficient) parameter, which seems to be a typo for B1.
dataframe = experiment.get('a1')

# %%
# Plotting the results
# Visualizing how the B1 (a3) parameter varies with the core diameter.
dataframe.plot_data(x="core_diameter")