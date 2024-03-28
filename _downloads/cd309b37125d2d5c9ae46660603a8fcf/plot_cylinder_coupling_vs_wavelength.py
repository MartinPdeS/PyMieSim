"""
Cylinder: Coupling vs Wavelength
================================

This example demonstrates how to compute and visualize the coupling efficiency as a function of wavelength for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment.detector import LPMode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim import measure
from PyMieSim.materials import BK7

# %%
# Defining the source
source_set = Gaussian(
    wavelength=np.linspace(950e-9, 1050e-9, 300),  # Wavelengths ranging from 950 nm to 1050 nm
    polarization_value=0,  # Linear polarization angle in radians
    polarization_type='linear',
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
# Here we look at cylinders with a set diameter, refractive index, and medium.
scatterer_set = Cylinder(
    diameter=np.linspace(100e-9, 8000e-9, 5),  # Diameters ranging from 100 nm to 8000 nm
    material=BK7,  # Material of the cylinder
    n_medium=1,  # Refractive index of the surrounding medium
    source_set=source_set
)

# %%
# Defining the detector
detector_set = LPMode(
    mode_number="LP11",  # Specifying the LP11 mode
    NA=[0.05, 0.01],  # Array of Numerical Apertures for the detector
    phi_offset=-180,  # Phi offset in degrees
    gamma_offset=0,  # Gamma offset in degrees
    polarization_filter=None,  # No polarization filter
    sampling=300  # Number of sampling points
)

# %%
# Setting up the experiment
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set,
    detector_set=detector_set
)

# %%
# Measuring the coupling efficiency
data = experiment.get(measure.coupling)

# %%
# Plotting the results
# Visualizing how the coupling efficiency varies with the wavelength.
figure = data.plot(
    x=source_set.wavelength,  # Wavelength as the x-axis
    std=scatterer_set.diameter  # Standard deviation with respect to cylinder diameter
)

# %%
# Displaying the plot
_ = figure.show()
