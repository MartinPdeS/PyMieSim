"""
Cylinder: Coupling vs Wavelength
================================

This example demonstrates how to compute and visualize the coupling efficiency as a function of wavelength for cylindrical scatterers using PyMieSim.
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

# %%
# Defining the source
source = Gaussian(
    wavelength=np.linspace(950, 1050, 300) * ureg.nanometer,  # Wavelengths ranging from 950 nm to 1050 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1e-3 * ureg.watt,  # 1 milliureg.watt
    NA=0.2 * ureg.AU  # Numerical Aperture
)

# %%
# Defining the scatterer distribution
# Here we look at cylinders with a set diameter, refractive index, and medium.
scatterer = Cylinder(
    diameter=np.linspace(100, 8000, 5) * ureg.nanometer,  # Diameters ranging from 100 nm to 8000 nm
    property=Material.BK7,  # Material of the cylinder
    medium_property=1 * ureg.RIU,  # Refractive index of the surrounding medium
    source=source
)

# %%
# Defining the detector
detector = CoherentMode(
    mode_number="LP11",  # Specifying the LP11 mode
    NA=[0.05, 0.01] * ureg.AU,  # Array of Numerical Apertures for the detector
    phi_offset=-180 * ureg.degree,  # Phi offset in ureg.degrees
    gamma_offset=0 * ureg.degree,  # Gamma offset in ureg.degrees
    polarization_filter=None,  # No polarization filter
    sampling=300 * ureg.AU,  # Number of sampling points
    rotation=0 * ureg.degree,  # Rotation of the mode field
)

# %%
# Setting up the experiment
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# %%
# Measuring the coupling efficiency
dataframe = experiment.get('coupling', scale_unit=True)

# %%
# Plotting the results
# Visualizing how the coupling efficiency varies with the wavelength.
dataframe.plot(x="source:wavelength", std='scatterer:diameter')
