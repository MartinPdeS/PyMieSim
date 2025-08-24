"""
CoreShell: B1 vs Core Diameter
==============================

This example demonstrates how to compute and visualize the B1 scattering parameter as a function of core diameter for CoreShell scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

# %%
# Defining the source to be employed.
# The source is always a plane wave in the LMT framework.
# The amplitude is set to one per default.
source = Gaussian(
    wavelength=800 * ureg.nanometer,  # 800 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU  # Numerical Aperture
)

# %%
# Defining the ranging parameters for the scatterer distribution
# Here we look at core/shell scatterers and use constant shell diameter
# with variable core diameter
scatterer = CoreShell(
    core_diameter=np.geomspace(100, 3000, 500) * ureg.nanometer,  # Geometrically spaced core diameters
    shell_thickness=800 * ureg.nanometer,  # Shell width of 800 nm
    core_property=1.6 * ureg.RIU,  # Refractive index of the core
    shell_property=Material.BK7,  # BK7 glass material for the shell
    medium_property=1 * ureg.RIU,  # Refractive index of the surrounding medium
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
dataframe.plot(x="scatterer:core_diameter")
