"""
PyMieSim vs PyMieScatt for Core-Shell Particles
===============================================
"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure

# PyMieScatt import
import PyMieScatt as pms

# Define parameters
wavelength = 600e-9  # Light source wavelength in meters
polarization_value = 0
polarization_type = 'linear'
optical_power = 1  # Power in watts
NA = 0.2  # Numerical aperture
medium_index = 1.2
core_index = 1.5
shell_index = 1.4
shell_width = 600e-9  # Shell width in meters
core_diameters = np.geomspace(10e-9, 500e-9, 400)  # Core diameters in meters

# Setup source
source = Gaussian(
    wavelength=wavelength,
    polarization=polarization_value,
    optical_power=optical_power,
    NA=NA
)

# Setup scatterer
scatterer = CoreShell(
    core_diameter=core_diameters,
    shell_width=shell_width,
    core_index=core_index,
    shell_index=shell_index,
    medium_index=medium_index,
    source=source
)

# Define experiment setup
experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=None  # No detector set specified
)

# Simulate using PyMieSim
sim_data = experiment.get(measure.Qsca, export_as_numpy=True).squeeze()

# Simulate using PyMieScatt
scatt_data = np.array([
    pms.MieQCoreShell(
        mCore=core_index,
        mShell=shell_index,
        wavelength=wavelength,
        dCore=diameter,
        dShell=diameter + shell_width,
        nMedium=medium_index
    )[1]
    for diameter in core_diameters
])

# Plot results
plt.figure(figsize=(8, 4))
plt.plot(core_diameters * 1e6, sim_data, 'C1-', linewidth=3, label='PyMieSim')
plt.plot(core_diameters * 1e6, scatt_data, 'k--', linewidth=1, label='PyMieScatt')
plt.xlabel('Core Diameter (μm)')
plt.ylabel('Scattering Efficiency')
plt.title('Scattering Efficiency Comparison for Core-Shell Particles')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
