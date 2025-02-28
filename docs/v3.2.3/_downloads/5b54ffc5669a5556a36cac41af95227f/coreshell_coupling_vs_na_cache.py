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
from PyOptik import Material
from PyMieSim.units import micrometer, nanometer, degree, watt, AU, RIU

# %%
# Defining the source
source = Gaussian(
    wavelength=1.2 * micrometer,  # 1200 nm
    polarization=90 * degree,  # Polarization angle in degrees
    optical_power=1e-3 * watt,  # 1 milliwatt
    NA=0.2 * AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
scatterer = CoreShell(
    core_diameter=1500 * nanometer,  # Core diameters from 100 nm to 600 nm
    shell_thickness=800 * nanometer,  # Shell width of 800 nm
    core_property=Material.silver,  # Core material
    shell_property=Material.BK7,  # Shell material
    medium_property=1 * RIU,  # Surrounding medium's refractive index
    source=source
)

# %%
# Defining the detector
detector = Photodiode(
    NA=[0.3] * AU,  # Numerical Apertures for the detector
    cache_NA=numpy.linspace(0.0, 0.2, 200) * AU,
    phi_offset=-180.0 * degree,  # Phi offset in degrees
    gamma_offset=0.0 * degree,  # Gamma offset in degrees
    sampling=4000 * AU,  # Number of sampling points
    polarization_filter=1 * degree
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the coupling efficiency
dataframe = experiment.get('coupling')

# %%
# Plotting the results
# Visualizing how the coupling efficiency varies with the core diameter.
dataframe.plot(x="detector:cache_NA")
