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
from PyOptik import UsualMaterial
from PyMieSim.experiment import measure

# %%
# Defining the source
source = Gaussian(
    wavelength=1.2e-6,  # 1200 nm
    polarization=90,  # Polarization angle in degrees
    optical_power=1e-3,  # 1 milliwatt
    NA=0.2  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = Cylinder(
    diameter=2e-6,  # 2000 nm
    material=UsualMaterial.BK7,  # Material of the cylinder
    medium_index=1,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Defining the detector
detector = Photodiode(
    NA=[0.5, 0.3, 0.1, 0.05],  # Array of Numerical Apertures for the detector
    phi_offset=np.linspace(-180, 180, 400),  # Angular displacement from -180 to 180 degrees
    gamma_offset=0,  # Gamma offset in degrees
    sampling=400,  # Number of sampling points
    polarization_filter=None  # No polarization filter
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
# Visualizing how the coupling efficiency varies with angular displacement.
figure = data.plot(
    x=detector.phi_offset,  # Angular displacement as the x-axis
    y_scale='log',  # Logarithmic scale for the y-axis
    normalize=True  # Normalizing the results
)

# %%
# Displaying the plot
_ = figure.show()
