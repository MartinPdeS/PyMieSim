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
from PyMieSim.materials import BK7
from PyMieSim import measure

# %%
# Defining the source
source_set = Gaussian(
    wavelength=1.2e-6,  # 1200 nm
    polarization_value=90,  # Polarization angle in degrees
    polarization_type='linear',
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer_set = Cylinder(
    diameter=2e-6,  # 2000 nm
    material=BK7,  # Material of the cylinder
    n_medium=1,  # Refractive index of the surrounding medium
    source_set=source_set
)

# %%
# Defining the detector
detector_set = Photodiode(
    NA=[0.5, 0.3, 0.1, 0.05],  # Array of Numerical Apertures for the detector
    phi_offset=np.linspace(-180, 180, 400),  # Angular displacement from -180 to 180 degrees
    gamma_offset=0,  # Gamma offset in degrees
    sampling=400,  # Number of sampling points
    polarization_filter=None  # No polarization filter
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
# Visualizing how the coupling efficiency varies with angular displacement.
figure = data.plot(
    x=detector_set.phi_offset,  # Angular displacement as the x-axis
    y_scale='log',  # Logarithmic scale for the y-axis
    normalize=True  # Normalizing the results
)

# %%
# Displaying the plot
_ = figure.show()
