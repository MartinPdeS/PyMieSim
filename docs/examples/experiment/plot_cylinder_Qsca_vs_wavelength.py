"""
Cylinder: Qsca vs Wavelength
============================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of wavelength for cylindrical scatterers using PyMieSim, considering cylinders with different diameters and refractive indices.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source
# Studying the scattering efficiency across a range of wavelengths.
source = Gaussian(
    wavelength=np.linspace(400, 1000, 150) * nanometer,  # Wavelengths ranging from 400 nm to 1000 nm
    polarization=0 * degree,  # Linear polarization angle in radians
    optical_power=1e-3 * watt,  # 1 milliwatt
    NA=0.2 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
# Considering cylinders with specific diameters and refractive indices.
scatterer = Cylinder(
    diameter=[200, 150] * nanometer,  # Array of diameters: 200 nm, 150 nm, 100 nm
    index=[2, 3, 4] * RIU,  # Array of refractive indices: 2, 3, 4
    medium_index=[1] * RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the scattering efficiency (Qsca)
# Averaging the data across the different indices to simplify visualization.
dataframe = experiment.get('Qsca')

# %%
# Plotting the results
# Visualizing how the Qsca varies with wavelength for the given cylinder configurations.
dataframe.plot_data(x='wavelength', std='index')
