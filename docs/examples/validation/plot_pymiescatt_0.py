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
from PyMieSim.units import degree, watt, AU, RIU, nanometer
from PyMieSim.utils import get_pymiescatt_sphere_dataframe
from MPSPlots.styles import mps

# PyMieScatt import
import PyMieScatt as pms

# Define parameters
wavelength = 632.8 * nanometer  # Wavelength of the light source in meters
polarization_value = 0 * degree
optical_power = 1e-3 * watt  # Power in watts
NA = 0.2 * AU  # Numerical aperture
medium_index = 1.21 * RIU
index = 1.4 * RIU
diameters = np.geomspace(10, 1_000, 50) * nanometer  # Geometric space for diameters

# Setup source
source = Gaussian(
    wavelength=wavelength,
    polarization=polarization_value,
    optical_power=optical_power,
    NA=NA
)

# Setup scatterer
scatterer = Sphere(
    diameter=diameters,
    index=index,
    medium_index=medium_index,
    source=source
)

# Define experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# Simulate using PyMieSim
pymiesim_dataframe = experiment.get('Qsca', 'Qext', 'Qabs', 'g', 'Qpr', 'Qback').reset_index('diameter')

pymiescatt_dataframe = get_pymiescatt_sphere_dataframe(
    wavelengths=wavelength,
    indexes=index,
    diameters=diameters,
    medium_indexes=medium_index
).reset_index('diameter')

# Plot results
with plt.style.context(mps):
    figure, ax = plt.subplots(1, 1)
    pymiescatt_dataframe.plot(x='diameter', ax=ax, linewidth=3)
    pymiesim_dataframe.plot(x='diameter', ax=ax, color='black', linestyle='--', linewidth=1.5)
    ax.set(
        xlabel='Diameter [μm]',
        ylabel='Scattering Efficiency',
        title='Scattering Efficiency Comparison for Sphere Particles'
    )
    plt.legend()
    plt.show()

