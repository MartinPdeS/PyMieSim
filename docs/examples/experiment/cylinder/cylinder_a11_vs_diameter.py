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
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source
source = Gaussian(
    wavelength=400 * nanometer,  # 400 nm
    polarization=90 * degree,  # Polarization angle in degrees
    optical_power=1e-3 * watt,  # 1 milliwatt
    NA=0.2 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = Cylinder(
    diameter=np.linspace(100, 10000, 800) * nanometer,  # Diameters ranging from 100 nm to 10000 nm
    property=1.4 * RIU,  # Refractive index of the cylinder
    medium_property=[1, 1.1] * RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the A1 scattering coefficient
# Note: The original request was for "a21"; assuming it meant A1, as "a21" might be a typo.
dataframe = experiment.get('a21')

# %%
# Plotting the results
# Visualizing how the A1 scattering coefficient varies with the cylinder diameter.
dataframe.plot_data(x='scatterer:diameter')