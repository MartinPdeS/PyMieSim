"""
InfiniteCylinder Scatterer Bohren-Huffman figure 8.8
====================================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
from PyMieSim.units import ureg

# PyMieSim imports
from PyMieSim.directories import validation_data_path
from PyMieSim.experiment.scatterer import InfiniteCylinderSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup

# Load theoretical data
theoretical_data = np.genfromtxt(
    f"{validation_data_path}/bohren_huffman/figure_88.csv", delimiter=","
)

# Define parameters
wavelength = 632.8 * ureg.nanometer  # Wavelength of the source in meters
polarization_set = PolarizationSet(angles=[0, 90] * ureg.degree)  # Create polarization set
optical_power = 1e-3 * ureg.watt  # Optical power in watts
numerical_aperture = 0.2 * ureg.AU  # Numerical aperture
diameters = np.geomspace(10, 6000, 800) * ureg.nanometer  # Diameters from 10 nm to 6 μm
index = 1.55 * ureg.RIU  # Refractive index of the cylinder
medium_index = 1.335 * ureg.RIU  # Refractive index of the medium

# Calculate the volume of the cylinders
volumes = np.pi * (diameters / 2) ** 2

# Configure the Gaussian source
source = GaussianSet(
    wavelength=[632.8] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

# Setup cylindrical scatterers
scatterer = InfiniteCylinderSet(
    diameter=diameters,
    refractive_index=[1.55] * ureg.RIU,
    medium_refractive_index=[1.335] * ureg.RIU,
    source=source
)

# Create experimental setup
experiment = Setup(scatterer=scatterer, source=source)

# Compute PyMieSim scattering cross section data
csca_data = experiment.get("Csca", as_numpy=True)

normalized_csca = (
    csca_data / volumes.to_base_units() * 1e-4 / 100
)  # Normalize the data as per specific needs

# Plotting the results
plt.figure(figsize=(8, 4))
plt.plot(
    diameters * 1e6,
    normalized_csca[0],
    "C0-",
    linewidth=3,
    label="PyMieSim Polarization: 0",
)
plt.plot(
    diameters * 1e6,
    normalized_csca[1],
    "C1-",
    linewidth=3,
    label="PyMieSim Polarization: 90",
)
plt.plot(
    diameters * 1e6,
    theoretical_data[0],
    "k--",
    linewidth=1,
    label="Theoretical BH 8.8 Polarization: 0",
)
plt.plot(
    diameters * 1e6,
    theoretical_data[1],
    "k--",
    linewidth=1,
    label="Theoretical BH 8.8 Polarization: 90",
)

plt.xlabel("Diameter (μm)")
plt.ylabel("Normalized Scattering Cross Section")
plt.title("Comparison of Scattering Cross Sections for InfiniteCylinders")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
