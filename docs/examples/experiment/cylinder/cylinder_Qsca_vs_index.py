"""
Cylinder: Qsca vs Index
=======================

This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a function of refractive index for cylindrical scatterers using PyMieSim, considering multiple wavelengths.
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
    wavelength=[500, 1000, 1500]
    * ureg.nanometer,  # Array of wavelengths: 500 nm, 1000 nm, 1500 nm
    polarization=30 * ureg.degree,  # Polarization angle in ureg.degrees
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU,  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = Cylinder(
    diameter=800 * ureg.nanometer,  # Fixed diameter of 800 nm
    property=np.linspace(1.3, 1.9, 1500)
    * ureg.RIU,  # Refractive index ranging from 1.3 to 1.9
    medium_property=1 * ureg.RIU,  # Refractive index of the surrounding medium
    source=source,
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the scattering efficiency (Qsca)
dataframe = experiment.get("Qsca", "Qext")

# %%
# Plotting the results
# Visualizing how the Qsca varies with the refractive index of the cylinder.
dataframe.plot(x="scatterer:property")
