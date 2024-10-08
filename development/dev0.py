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
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source
source = Gaussian(
    wavelength=numpy.geomspace(300, 1200, 400) * nanometer,  # Array of wavelengths: 800 nm, 900 nm, 1000 nm
    polarization=0 * degree,  # Linear polarization angle in radians
    optical_power=1e-3 * watt,  # 1 milliwatt
    NA=0.2 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = CoreShell(
    core_diameter=600 * nanometer,  # Core diameters from 100 nm to 600 nm
    shell_width=numpy.linspace(10, 400, 20) * nanometer,  # Shell width of 800 nm
    core_property=[Material.silver, Material.gold],  # Core material
    shell_property=Material.BK7,  # Shell material
    medium_property=1 * RIU,  # Surrounding medium's refractive index
    source=source
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the backscattering efficiency (Qback)
# For demonstrating the measurement of Qsca, a separate call to `experiment.get()` with `measure.Qsca` is needed.
dataframe = experiment.get('Qsca')

# %%
# Plotting the results
# Visualizing how the backscattering efficiency varies with the core diameter.
dataframe.plot_data(x="wavelength", std='shell_width')