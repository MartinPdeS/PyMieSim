"""
InfiniteCylinder Scatterer Bohren-Huffman figure 8.8
====================================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
from PyMieSim import (
    Measure,
    ureg,
    InfiniteCylinderSet,
    GaussianSet,
    PolarizationSet,
    Experiment,
)

# PyMieSim imports
from PyMieSim.directories import validation_data_path

# Load theoretical data
theoretical_data = np.genfromtxt(
    f"{validation_data_path}/bohren_huffman/figure_88.csv", delimiter=","
)

# Define parameters
medium_index = 1.335
# The packaged reference curves use 632.8 nm as the wavelength in the medium.
vacuum_wavelength = medium_index * (632.8 * ureg.nanometer)
polarization_set = PolarizationSet(angles=[0, 90] * ureg.degree)  # Create polarization set
optical_power = 1e-3 * ureg.watt  # Optical power in watts
numerical_aperture = 0.2  # Numerical aperture
diameters = np.geomspace(10, 6000, 800) * ureg.nanometer  # Diameters from 10 nm to 6 μm
index = 1.55  # Refractive index of the cylinder

# Cross-sectional area is cylinder volume per unit axial length.
geometric_areas = np.pi * (diameters / 2) ** 2

# Configure the Gaussian source
source = GaussianSet(
    wavelength=vacuum_wavelength,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

# Setup cylindrical scatterers
scatterer = InfiniteCylinderSet(
    diameter=diameters,
    material=[1.55],
    medium=[medium_index],
)

# Create experimental setup
experiment = Experiment(scatterer_set=scatterer, source_set=source)

# Compute PyMieSim scattering cross section data
csca_data = experiment.run(Measure.CSCA).quantity

normalized_csca = (
    (csca_data / geometric_areas).to("1 / micrometer").magnitude
)  # Normalize the data as per specific needs

# Plotting the results
plt.figure(figsize=(8, 4))
plt.plot(
    diameters.to("micrometer").magnitude,
    normalized_csca[0],
    "C0-",
    linewidth=3,
    label="PyMieSim Polarization: 0",
)
plt.plot(
    diameters.to("micrometer").magnitude,
    normalized_csca[1],
    "C1-",
    linewidth=3,
    label="PyMieSim Polarization: 90",
)
plt.plot(
    diameters.to("micrometer").magnitude,
    theoretical_data[0],
    "k--",
    linewidth=1,
    label="Theoretical BH 8.8 Polarization: 0",
)
plt.plot(
    diameters.to("micrometer").magnitude,
    theoretical_data[1],
    "k--",
    linewidth=1,
    label="Theoretical BH 8.8 Polarization: 90",
)

plt.xlabel("Diameter (μm)")
plt.ylabel("Scattering width / geometric area (μm⁻¹)")
plt.title("Comparison of Scattering Cross Sections for InfiniteCylinders")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
