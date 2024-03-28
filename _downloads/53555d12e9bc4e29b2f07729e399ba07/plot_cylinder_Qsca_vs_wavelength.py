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
from PyMieSim import measure

# %%
# Defining the source
# Studying the scattering efficiency across a range of wavelengths.
source_set = Gaussian(
    wavelength=np.linspace(400e-9, 1000e-9, 500),  # Wavelengths ranging from 400 nm to 1000 nm
    polarization_value=0,  # Linear polarization angle in radians
    polarization_type='linear',
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
# Considering cylinders with specific diameters and refractive indices.
scatterer_set = Cylinder(
    diameter=[200e-9, 150e-9, 100e-9],  # Array of diameters: 200 nm, 150 nm, 100 nm
    index=[2, 3, 4],  # Array of refractive indices: 2, 3, 4
    n_medium=1,  # Refractive index of the surrounding medium
    source_set=source_set
)

# %%
# Setting up the experiment
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the scattering efficiency (Qsca)
# Averaging the data across the different indices to simplify visualization.
data = experiment.get(measure.Qsca)
data = data.mean(scatterer_set.index)

# %%
# Plotting the results
# Visualizing how the Qsca varies with wavelength for the given cylinder configurations.
figure = data.plot(
    x=source_set.wavelength,  # Wavelength as the x-axis
    y_scale='linear'  # Linear scale for the y-axis
)

# %%
# Displaying the plot
_ = figure.show()
