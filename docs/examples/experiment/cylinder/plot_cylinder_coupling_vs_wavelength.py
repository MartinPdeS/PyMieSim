"""
Cylinder: Coupling vs Wavelength
================================

This example demonstrates how to compute and visualize the coupling efficiency as a function of wavelength for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure
from PyOptik import UsualMaterial

# %%
# Defining the source
source = Gaussian(
    wavelength=np.linspace(950e-9, 1050e-9, 300),  # Wavelengths ranging from 950 nm to 1050 nm
    polarization=0,  # Linear polarization angle in radians
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
# Here we look at cylinders with a set diameter, refractive index, and medium.
scatterer = Cylinder(
    diameter=np.linspace(100e-9, 8000e-9, 5),  # Diameters ranging from 100 nm to 8000 nm
    material=UsualMaterial.BK7,  # Material of the cylinder
    medium_index=1,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Defining the detector
detector = CoherentMode(
    mode_number="LP11",  # Specifying the LP11 mode
    NA=[0.05, 0.01],  # Array of Numerical Apertures for the detector
    phi_offset=-180,  # Phi offset in degrees
    gamma_offset=0,  # Gamma offset in degrees
    polarization_filter=None,  # No polarization filter
    sampling=300,  # Number of sampling points
    rotation=0,  # Rotation of the mode field
)

# %%
# Setting up the experiment
experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

# %%
# Measuring the coupling efficiency
data = experiment.get(measure.coupling)

# %%
# Plotting the results
# Visualizing how the coupling efficiency varies with the wavelength.
figure = data.plot(
    x=source.wavelength,  # Wavelength as the x-axis
    std=scatterer.diameter  # Standard deviation with respect to cylinder diameter
)

# %%
# Displaying the plot
_ = figure.show()
