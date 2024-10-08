"""
Cylinder: Qsca vs Index
=======================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of refractive index for cylindrical scatterers using PyMieSim, considering multiple wavelengths.
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
    wavelength=[500, 1000, 1500] * nanometer,  # Array of wavelengths: 500 nm, 1000 nm, 1500 nm
    polarization=30 * degree,  # Polarization angle in degrees
    optical_power=1e-3 * watt,  # 1 milliwatt
    NA=0.2 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = Cylinder(
    diameter=800 * nanometer,  # Fixed diameter of 800 nm
    property=np.linspace(1.3, 1.9, 1500) * RIU,  # Refractive index ranging from 1.3 to 1.9
    medium_property=1 * RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the scattering efficiency (Qsca)
dataframe = experiment.get('Qsca', 'Qext')

# %%
# Plotting the results
# Visualizing how the Qsca varies with the refractive index of the cylinder.
dataframe.plot_data(x="property")
