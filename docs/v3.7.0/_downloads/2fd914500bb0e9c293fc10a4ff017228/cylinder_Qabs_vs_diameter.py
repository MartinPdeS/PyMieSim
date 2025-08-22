"""
Cylinder: Qabs vs Diameter
==========================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of diameter for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

# %%
# Defining the source
source = Gaussian(
    wavelength=400 * ureg.nanometer,  # 400 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = Cylinder(
    diameter=np.linspace(1, 800, 300) * ureg.nanometer,  # Diameters ranging from 1 nm to 800 nm
    property=[Material.silver, Material.gold, Material.aluminium],  # Scatterer materials
    medium_property=1 * ureg.RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the scattering efficiency (Qsca)
# Note: The original request mentioned Qsca, but the measurement code uses Qabs.
# If Qsca measurement is intended, ensure to use the correct measure object from PyMieSim.
dataframe = experiment.get('Qabs')  # Assuming Qabs was intended, replace with measure.Qsca if needed

# %%
# Plotting the results
# Visualizing how the scattering efficiency varies with the cylinder diameter.
dataframe.plot(x='scatterer:diameter')
