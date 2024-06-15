"""
Comparison of Scattering Efficiency Using PyMieSim vs PyMieScatt
================================================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure

# PyMieScatt import
import PyMieScatt as pms

# Define parameters
wavelength = 632.8e-9  # Wavelength of the source in meters
index = 1.4  # Refractive index of the sphere
medium_index = 1.0  # Refractive index of the medium
optical_power = 1  # Power of the light source in watts
NA = 0.2  # Numerical aperture
diameters = np.geomspace(10e-9, 6e-6, 800)  # Diameters from 10 nm to 6 μm

# Configure the Gaussian source
source = Gaussian(
    wavelength=wavelength,
    polarization=0,
    optical_power=optical_power,
    NA=NA
)

# Setup spherical scatterer
scatterer = Sphere(
    diameter=diameters,
    index=index,
    medium_index=medium_index,
    source=source
)

# Create experimental setup
experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=None  # No detector configuration
)

# Compute PyMieSim scattering efficiency data
sim_data = experiment.get(measure.Qsca, export_as_numpy=True).squeeze()

# Compute PyMieScatt scattering efficiency data
scatt_data = np.array([
    pms.MieQ(
        m=index,
        wavelength=wavelength,
        diameter=d
    )[1] for d in diameters
]).squeeze()

# Plotting the results
plt.figure(figsize=(8, 4))
plt.plot(diameters * 1e6, sim_data, 'C1-', linewidth=3, label='PyMieSim')
plt.plot(diameters * 1e6, scatt_data, 'k--', linewidth=1, label='PyMieScatt')
plt.xlabel('Diameter (μm)')
plt.ylabel('Scattering Efficiency')
plt.title('Scattering Efficiency Comparison: Sphere')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
