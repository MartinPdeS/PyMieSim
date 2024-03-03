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
from PyMieSim import measure
from PyMieSim.materials import BK7

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
    diameter=np.linspace(100e-9, 3000e-9, 200),  # Diameters ranging from 100 nm to 3000 nm
    material=BK7,  # Material of the cylinder
    n_medium=1.0,  # Refractive index of the surrounding medium
    source_set=source_set
)

# %%
# Defining the detector
detector_set = Photodiode(
    NA=[0.1, 0.05],  # Numerical Apertures for the detector
    phi_offset=-180.0,  # Phi offset in degrees
    gamma_offset=0.0,  # Gamma offset in degrees
    sampling=600,  # Number of sampling points
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
# Visualizing how the coupling efficiency varies with the cylinder diameter.
figure = data.plot(
    x=scatterer_set.diameter,  # Cylinder diameter as the x-axis
    y_scale='linear',  # Linear scale for the y-axis
    normalize=True  # Normalizing the results
)

# %%
# Displaying the plot
_ = figure.show()
