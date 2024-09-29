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
from PyMieSim.units import degree, watt, AU, RIU, nanometer

# PyMieScatt import
import PyMieScatt as pms

# Define parameters
wavelength = 600 * nanometer  # Light source wavelength in meters
polarization = 0 * degree
optical_power = 1 * watt  # Power in watts
NA = 0.2 * AU  # Numerical aperture
medium_index = 1.2 * RIU
core_index = 1.5 * RIU
shell_index = 1.4 * RIU
shell_width = 600 * nanometer  # Shell width in meters
core_diameters = np.geomspace(10, 500, 400) * nanometer  # Core diameters in meters

# Setup source
source = Gaussian(
    wavelength=wavelength,
    polarization=polarization,
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
sim_data = experiment.get(measure.Qsca, export_as='numpy').squeeze()

# Simulate using PyMieScatt
scatt_data = np.array([
    pms.MieQCoreShell(
        mCore=core_index.to_base_units().magnitude,
        mShell=shell_index.to_base_units().magnitude,
        wavelength=wavelength.to_base_units().magnitude,
        dCore=diameter,
        dShell=diameter + shell_width.to_base_units().magnitude,
        nMedium=medium_index.to_base_units().magnitude
    )[1]
    for diameter in core_diameters.to_base_units().magnitude
])

# Plot results
plt.figure(figsize=(8, 4))
plt.plot(core_diameters, sim_data, 'C1-', linewidth=3, label='PyMieSim')
plt.plot(core_diameters, scatt_data, 'k--', linewidth=1, label='PyMieScatt')
plt.xlabel('Core Diameter [μm]')
plt.ylabel('Scattering Efficiency')
plt.title('Scattering Efficiency Comparison for Core-Shell Particles')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
