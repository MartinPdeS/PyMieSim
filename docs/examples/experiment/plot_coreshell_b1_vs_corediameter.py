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
from PyMieSim.materials import BK7
from PyMieSim import measure

# %%
# Defining the source to be employed.
# The source is always a plane wave in the LMT framework.
# The amplitude is set to one per default.
source_set = Gaussian(
    wavelength=800e-9,  # 800 nm
    polarization_value=0,  # Linear polarization angle in radians
    polarization_type='linear',
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the ranging parameters for the scatterer distribution
# Here we look at core/shell scatterers and use constant shell diameter
# with variable core diameter
scatterer_set = CoreShell(
    core_diameter=np.geomspace(100e-9, 3000e-9, 5000),  # Geometrically spaced core diameters
    shell_width=800e-9,  # Shell width of 800 nm
    core_index=1.6,  # Refractive index of the core
    shell_material=BK7,  # BK7 glass material for the shell
    n_medium=1,  # Refractive index of the surrounding medium
    source_set=source_set
)

# %%
# Defining the experiment setup
# Integrating the defined source and scatterers into a single experimental setup.
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the B1 scattering parameter
data = experiment.get(measure.b1)

# %%
# Plotting the results
# Visualizing how the B1 parameter varies with the core diameter.
figure = data.plot(
    x=scatterer_set.core_diameter,  # Core diameter as the x-axis
    y_scale='linear'  # Linear scale for the y-axis
)

# %%
# Displaying the plot
_ = figure.show()
