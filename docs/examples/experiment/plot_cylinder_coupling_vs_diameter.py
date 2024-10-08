"""
Cylinder: Coupling vs Diameter
==============================

This example demonstrates how to compute and visualize the coupling efficiency as a function of diameter for cylindrical scatterers using PyMieSim.
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
    NA=0.2 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = Cylinder(
    diameter=np.linspace(100, 3000, 200) * nanometer,  # Diameters ranging from 100 nm to 3000 nm
    property=[Material.BK7],  # Material of the cylinder
    medium_property=1.0 * RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Defining the detector
detector = Photodiode(
    NA=[0.1, 0.05] * AU,  # Numerical Apertures for the detector
    phi_offset=-180.0 * degree,  # Phi offset in degrees
    gamma_offset=0.0 * degree,  # Gamma offset in degrees
    sampling=600 * AU,  # Number of sampling points
    polarization_filter=None  # No polarization filter
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the coupling efficiency
dataframe = experiment.get('coupling')

# %%
# Plotting the results
# Visualizing how the coupling efficiency varies with the cylinder diameter.
# dataframe.plot_data(x="diameter")