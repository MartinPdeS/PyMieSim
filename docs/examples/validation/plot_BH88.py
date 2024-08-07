"""
Comparison of PyMieSim and Theoretical Bohren-Huffman Data for Cylinder Scattering
==================================================================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.directories import validation_data_path
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.experiment import measure

# Load theoretical data
theoretical_data = np.genfromtxt(f"{validation_data_path}/Figure88BH.csv", delimiter=',')

# Define parameters
wavelength = 632.8e-9  # Wavelength of the source in meters
polarization_values = [0, 90]  # Polarization values in degrees
optical_power = 1e-3  # Optical power in watts
NA = 0.2  # Numerical aperture
diameters = np.geomspace(10e-9, 6e-6, 800)  # Diameters from 10 nm to 6 μm
index = 1.55  # Refractive index of the cylinder
medium_index = 1.335  # Refractive index of the medium

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
    index=index,
    medium_index=medium_index,
    source=source
)

# Create experimental setup
experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=None
)

# Compute PyMieSim scattering cross section data
csca_data = experiment.get(measure.Csca, export_as_numpy=True).squeeze()
normalized_csca = csca_data / volumes * 1e-4 / 100  # Normalize the data as per specific needs

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

# Verify data accuracy
assert np.allclose(normalized_csca, theoretical_data, atol=1e-9), 'Error: mismatch on BH 8.8 calculation occurring'
