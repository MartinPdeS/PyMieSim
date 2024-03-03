"""
CoreShell: Qback vs Core Diameter
=======================================

This example demonstrates how to compute and visualize the backscattering efficiency (Qback)
as functions of core diameter for CoreShell scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.materials import BK7, Silver
from PyMieSim import measure

# %%
# Defining the source
source_set = Gaussian(
    wavelength=[800e-9, 900e-9, 1000e-9],  # Array of wavelengths: 800 nm, 900 nm, 1000 nm
    polarization_value=0,  # Linear polarization angle in radians
    polarization_type='linear',
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer_set = CoreShell(
    core_diameter=numpy.geomspace(100e-9, 600e-9, 400),  # Core diameters from 100 nm to 600 nm
    shell_width=800e-9,  # Shell width of 800 nm
    core_material=Silver,  # Core material
    shell_material=BK7,  # Shell material
    n_medium=1,  # Surrounding medium's refractive index
    source_set=source_set
)

# %%
# Setting up the experiment
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the backscattering efficiency (Qback)
# For demonstrating the measurement of Qsca, a separate call to `experiment.get()` with `measure.Qsca` is needed.
data = experiment.get(measure.Qback)

# %%
# Plotting the results
# Visualizing how the backscattering efficiency varies with the core diameter.
figure = data.plot(
    x=scatterer_set.core_diameter,  # Core diameter as the x-axis
    y_scale='log'  # Logarithmic scale for the y-axis
)

# %%
# Displaying the plot
_ = figure.show()
