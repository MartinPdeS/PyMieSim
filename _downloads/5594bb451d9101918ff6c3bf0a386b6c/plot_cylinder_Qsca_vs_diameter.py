"""
Cylinder: Qsca vs Diameter
==========================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim, considering multiple wavelengths.
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
source_set = Gaussian(
    wavelength=[500e-9, 1000e-9, 1500e-9],  # Array of wavelengths: 500 nm, 1000 nm, 1500 nm
    polarization_value=30,  # Polarization angle in degrees
    polarization_type='linear',
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer_set = Cylinder(
    diameter=np.geomspace(6.36e-9, 10000e-9, 1000),  # Diameters ranging from ~6.36 nm to 10000 nm
    index=[1.4],  # Refractive index of the cylinder
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
data = experiment.get(measure.Qsca)

# %%
# Plotting the results
# Visualizing how the Qsca varies with the cylinder diameter.
figure = data.plot(
    x=scatterer_set.diameter,  # Cylinder diameter as the x-axis
    y_scale='linear'  # Linear scale for the y-axis
)

# %%
# Displaying the plot
_ = figure.show()
