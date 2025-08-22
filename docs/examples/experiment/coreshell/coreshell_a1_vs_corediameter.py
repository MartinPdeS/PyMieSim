"""
CoreShell: An vs Core Diameter
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
# Defining the source
# In the LMT framework, the source is always considered a plane wave with a default amplitude of one.
source = Gaussian(
    wavelength=800 * ureg.nanometer,  # 800 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
# Here, we explore core/shell scatterers with a constant shell diameter and variable core diameter.
scatterer = CoreShell(
    core_diameter=np.geomspace(100, 600, 10) * ureg.nanometer,  # Geometrically spaced core diameters
    shell_thickness=150 * ureg.nanometer,  # Shell width of 800 nm
    core_property=[1.4] * ureg.RIU,  # Refractive index of the core
    shell_property=[Material.BK7],  # BK7 glass material for the shell
    medium_property=1 * ureg.RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Defining the experiment setup
# Integrating the defined source and scatterers into a single experimental setup.
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the B1 scattering parameter
# Here, we're interested in the a3 (first magnetic coefficient) parameter, which seems to be a typo for B1.
dataframe = experiment.get('a1', 'a2', 'a3',  'Qsca')

# %%
# Plotting the results
# Visualizing how the B1 (a3) parameter varies with the core diameter.
dataframe.plot(x="scatterer:core_diameter")
