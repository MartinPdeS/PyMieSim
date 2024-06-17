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
from PyMieSim.experiment import measure
from PyOptik import UsualMaterial

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
scatterer = CoreShell(
    core_diameter=numpy.geomspace(100e-9, 600e-9, 400),  # Core diameters from 100 nm to 600 nm
    shell_width=800e-9,  # Shell width of 800 nm
    core_material=UsualMaterial.Silver,  # Core material
    shell_material=UsualMaterial.BK7,  # Shell material
    medium_index=1,  # Surrounding medium's refractive index
    source=source
)

# %%
# Defining the detector
detector = Photodiode(
    NA=[0.1, 0.05],  # Numerical Apertures for the detector
    phi_offset=-180.0,  # Phi offset in degrees
    gamma_offset=0.0,  # Gamma offset in degrees
    sampling=600,  # Number of sampling points
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
# Visualizing how the coupling efficiency varies with the core diameter.
figure = data.plot(
    x=scatterer.core_diameter,  # Core diameter as the x-axis
    y_scale='linear',  # Linear scale for the y-axis
    normalize=True,  # Normalizing the results
)

# %%
# Displaying the plot
_ = figure.show()
