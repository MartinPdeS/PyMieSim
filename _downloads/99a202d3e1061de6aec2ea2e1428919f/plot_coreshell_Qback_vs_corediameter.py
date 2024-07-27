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
from PyOptik import UsualMaterial
from PyMieSim.experiment import measure

# %%
# Defining the source
source = Gaussian(
    wavelength=[800e-9, 900e-9, 1000e-9],  # Array of wavelengths: 800 nm, 900 nm, 1000 nm
    polarization=0,  # Linear polarization angle in radians
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = CoreShell(
    core_diameter=numpy.geomspace(100e-9, 600e-9, 400),  # Core diameters from 100 nm to 600 nm
    shell_width=800e-9,  # Shell width of 800 nm
    core_material=UsualMaterial.Silver,  # Core material
    shell_material=UsualMaterial.BK7,  # Shell material
    medium_index=1,  # Surrounding medium's refractive index
    source=source
)

# %%
# Setting up the experiment
experiment = Setup(
    scatterer=scatterer,
    source=source
)

# %%
# Measuring the backscattering efficiency (Qback)
# For demonstrating the measurement of Qsca, a separate call to `experiment.get()` with `measure.Qsca` is needed.
data = experiment.get(measure.Qback)

# %%
# Plotting the results
# Visualizing how the backscattering efficiency varies with the core diameter.
figure = data.plot(
    x=scatterer.core_diameter,  # Core diameter as the x-axis
    y_scale='log'  # Logarithmic scale for the y-axis
)

# %%
# Displaying the plot
_ = figure.show()
