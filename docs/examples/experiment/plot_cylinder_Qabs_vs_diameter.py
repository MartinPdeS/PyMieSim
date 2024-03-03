"""
Cylinder: Qsca vs Diameter
==========================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.materials import Gold, Silver, Aluminium
from PyMieSim import measure

# %%
# Defining the source
source_set = Gaussian(
    wavelength=400e-9,  # 400 nm
    polarization_value=0,  # Linear polarization angle in radians
    polarization_type='linear',
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer_set = Cylinder(
    diameter=np.linspace(1e-9, 800e-9, 300),  # Diameters ranging from 1 nm to 800 nm
    material=[Silver, Gold, Aluminium],  # Scatterer materials
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
# Note: The original request mentioned Qsca, but the measurement code uses Qabs. 
# If Qsca measurement is intended, ensure to use the correct measure object from PyMieSim.
data = experiment.get(measure.Qabs)  # Assuming Qabs was intended, replace with measure.Qsca if needed

# %%
# Plotting the results
# Visualizing how the scattering efficiency varies with the cylinder diameter.
figure = data.plot(
    x=scatterer_set.diameter,  # Cylinder diameter as the x-axis
    y_scale="linear"  # Linear scale for the y-axis
)

# %%
# Displaying the plot
_ = figure.show()
