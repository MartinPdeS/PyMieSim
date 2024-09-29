"""
PyMieSim vs PyMieScatt Comparison
=================================
"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure
from PyMieSim.units import degree, watt, AU, RIU, nanometer

# PyMieScatt import
import PyMieScatt as pms

# Define parameters
wavelength = 632.8 * nanometer  # Wavelength of the light source in meters
polarization_value = 0 * degree
optical_power = 1e-3 * watt  # Power in watts
NA = 0.2 * AU  # Numerical aperture
medium_index = 1.21 * RIU
refractive_index = 1.4 * RIU
diameter_range = np.geomspace(10, 1_000, 50) * nanometer  # Geometric space for diameters

# Setup source
source = Gaussian(
    wavelength=wavelength,
    polarization=polarization_value,
    optical_power=optical_power,
    NA=NA
)

# Setup scatterer
scatterer = Sphere(
    diameter=diameter_range,
    index=refractive_index,
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
    pms.MieQ(
        m=refractive_index.magnitude,
        diameter=d.to_base_units().magnitude,
        wavelength=wavelength.to_base_units().magnitude,
        nMedium=medium_index.magnitude)[1]
    for d in diameter_range
])

# Plot results
plt.figure(figsize=(8, 4))
plt.plot(diameter_range * 1e6, sim_data, 'C1-', linewidth=3, label='PyMieSim')
plt.plot(diameter_range * 1e6, scatt_data, 'k--', linewidth=1, label='PyMieScatt')
plt.xlabel('Diameter (μm)')
plt.ylabel('Scattering Efficiency')
plt.title('Comparison of Scattering Efficiency: PyMieSim vs PyMieScatt')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
