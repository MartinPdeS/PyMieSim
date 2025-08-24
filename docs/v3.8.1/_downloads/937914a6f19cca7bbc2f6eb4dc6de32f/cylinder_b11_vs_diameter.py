"""
Cylinder: B1 Scattering Coefficient
===================================

This example demonstrates how to compute and visualize the B1 scattering coefficient as a function of diameter for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

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
    diameter=np.linspace(100, 10000, 800) * ureg.nanometer,  # Diameters ranging from 100 nm to 10000 nm
    property=1.4 * ureg.RIU,  # Refractive index of the cylinder
    medium_property=1 * ureg.RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the B1 scattering coefficient
dataframe = experiment.get('b11')

# %%
# Plotting the results
# Visualizing how the B1 scattering coefficient varies with the cylinder diameter.
dataframe.plot(x="scatterer:diameter")
