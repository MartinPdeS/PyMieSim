"""
Cylinder Scatterer Bohren-Huffman figure 8.8
============================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.directories import validation_data_path
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import degree, watt, AU, RIU, nanometer

# Load theoretical data
theoretical_data = np.genfromtxt(f"{validation_data_path}/bohren_huffman/figure_88.csv", delimiter=',')

# Define parameters
wavelength = 632.8 * nanometer  # Wavelength of the source in meters
polarization_values = [0, 90] * degree  # Polarization values in degrees
optical_power = 1e-3 * watt  # Optical power in watts
NA = 0.2 * AU  # Numerical aperture
diameters = np.geomspace(10, 6000, 800) * nanometer  # Diameters from 10 nm to 6 μm
index = 1.55 * RIU  # Refractive index of the cylinder
medium_index = 1.335 * RIU  # Refractive index of the medium

# Calculate the volume of the cylinders
volumes = np.pi * (diameters / 2)**2

# Configure the Gaussian source
source = Gaussian(
    wavelength=wavelength,
    polarization=polarization_values,
    optical_power=optical_power,
    NA=NA
)

# Setup cylindrical scatterers
scatterer = Cylinder(
    diameter=diameters,
    property=index,
    medium_property=medium_index,
    source=source
)

# Create experimental setup
experiment = Setup(scatterer=scatterer, source=source)

# Compute PyMieSim scattering cross section data
csca_data = experiment.get('Csca', add_units=False).squeeze().values.reshape([-1, diameters.size])
normalized_csca = csca_data / volumes.to_base_units() * 1e-4 / 100  # Normalize the data as per specific needs

# Plotting the results
plt.figure(figsize=(8, 4))
plt.plot(diameters * 1e6, normalized_csca[0], 'C0-', linewidth=3, label='PyMieSim Polarization: 0')
plt.plot(diameters * 1e6, normalized_csca[1], 'C1-', linewidth=3, label='PyMieSim Polarization: 90')
plt.plot(diameters * 1e6, theoretical_data[0], 'k--', linewidth=1, label='Theoretical BH 8.8 Polarization: 0')
plt.plot(diameters * 1e6, theoretical_data[1], 'k--', linewidth=1, label='Theoretical BH 8.8 Polarization: 90')

plt.xlabel('Diameter (μm)')
plt.ylabel('Normalized Scattering Cross Section')
plt.title('Comparison of Scattering Cross Sections for Cylinders')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
