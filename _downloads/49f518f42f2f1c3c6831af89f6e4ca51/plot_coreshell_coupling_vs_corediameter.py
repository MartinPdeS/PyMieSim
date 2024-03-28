"""
CoreShell: Coupling vs Diameter
===============================

This example demonstrates how to compute and visualize the coupling efficiency as a function of core diameter for CoreShell scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim import measure
from PyMieSim.materials import BK7, Silver

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
scatterer_set = CoreShell(
    core_diameter=numpy.geomspace(100e-9, 600e-9, 400),  # Core diameters from 100 nm to 600 nm
    shell_width=800e-9,  # Shell width of 800 nm
    core_material=Silver,  # Core material
    shell_material=BK7,  # Shell material
    n_medium=1,  # Surrounding medium's refractive index
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
# Visualizing how the coupling efficiency varies with the core diameter.
figure = data.plot(
    x=scatterer_set.core_diameter,  # Core diameter as the x-axis
    y_scale='linear',  # Linear scale for the y-axis
    normalize=True,  # Normalizing the results
    add_box=True  # Adding a box around the plot for clarity
)

# %%
# Displaying the plot
_ = figure.show()
