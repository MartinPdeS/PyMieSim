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
from PyMieSim import measure

# %%
# Defining the source
source_set = Gaussian(
    wavelength=400e-9,  # 400 nm
    polarization_value=90,  # Polarization angle in degrees
    polarization_type='linear',
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer_set = Cylinder(
    diameter=np.linspace(100e-9, 10000e-9, 800),  # Diameters ranging from 100 nm to 10000 nm
    index=1.4,  # Refractive index of the cylinder
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
# Measuring the A1 scattering coefficient
# Note: The original request was for "a21"; assuming it meant A1, as "a21" might be a typo.
data = experiment.get(measure.a21)

# %%
# Plotting the results
# Visualizing how the A1 scattering coefficient varies with the cylinder diameter.
figure = data.plot(
    x=scatterer_set.diameter  # Cylinder diameter as the x-axis
)

# %%
# Displaying the plot
_ = figure.show()
