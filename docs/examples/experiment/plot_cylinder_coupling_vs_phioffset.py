"""
Cylinder: Goniometer
====================

This example demonstrates how to use a goniometer setup to measure and visualize the coupling efficiency as a function of angular displacement for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# %%
# Defining the source
source = Gaussian(
    wavelength=1200 * nanometer,  # 1200 nm
    polarization=90 * degree,  # Polarization angle in degrees
    optical_power=1e-3 * watt,  # 1 milliwatt
    NA=0.2 * AU # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = Cylinder(
    diameter=2000 * nanometer,  # 2000 nm
    property=Material.BK7,  # Material of the cylinder
    medium_property=1 * RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Defining the detector
detector = Photodiode(
    NA=[0.5, 0.3, 0.1, 0.05] * AU,  # Array of Numerical Apertures for the detector
    phi_offset=np.linspace(-180, 180, 200) * degree,  # Angular displacement from -180 to 180 degrees
    gamma_offset=0 * degree,  # Gamma offset in degrees
    sampling=400 * AU,  # Number of sampling points
    polarization_filter=None  # No polarization filter
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the coupling efficiency
dataframe = experiment.get('coupling', scale_unit=True)

# %%
# Plotting the results
# Visualizing how the coupling efficiency varies with angular displacement.
dataframe.plot_data(x="phi_offset")