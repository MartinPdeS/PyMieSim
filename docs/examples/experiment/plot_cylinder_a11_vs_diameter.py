"""
Cylinder: A1 Scattering Coefficient
===================================

This example demonstrates how to compute and visualize the A1 scattering coefficient as a function of diameter for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure

# %%
# Defining the source
source = Gaussian(
    wavelength=400e-9,  # 400 nm
    polarization=90,  # Polarization angle in degrees
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = Cylinder(
    diameter=np.linspace(100e-9, 10000e-9, 800),  # Diameters ranging from 100 nm to 10000 nm
    index=1.4,  # Refractive index of the cylinder
    medium_index=1,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Setting up the experiment
experiment = Setup(
    scatterer=scatterer,
    source=source
)

# %%
# Measuring the A1 scattering coefficient
# Note: The original request was for "a21"; assuming it meant A1, as "a21" might be a typo.
data = experiment.get(measure.a21)

# %%
# Plotting the results
# Visualizing how the A1 scattering coefficient varies with the cylinder diameter.
figure = data.plot(x=scatterer.diameter)

# %%
# Displaying the plot
_ = figure.show()
