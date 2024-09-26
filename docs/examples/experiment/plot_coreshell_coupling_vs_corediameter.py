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
    core_diameter=numpy.geomspace(100, 600, 400) * nanometer,  # Core diameters from 100 nm to 600 nm
    shell_width=800 * nanometer,  # Shell width of 800 nm
    core_material=Material.silver,  # Core material
    shell_material=Material.BK7,  # Shell material
    medium_index=1 * RIU,  # Surrounding medium's refractive index
    source=source
)

# %%
# Defining the detector
detector = Photodiode(
    NA=[0.1] * AU,  # Numerical Apertures for the detector
    phi_offset=-180.0 * degree,  # Phi offset in degrees
    gamma_offset=0.0 * degree,  # Gamma offset in degrees
    sampling=600 * AU,  # Number of sampling points
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
dataframe = experiment.get(measure.coupling)

# %%
# Plotting the results
# Visualizing how the coupling efficiency varies with the core diameter.
dataframe.plot_data(x="core_diameter")