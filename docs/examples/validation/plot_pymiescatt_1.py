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
from PyMieSim.units import degree, watt, AU, RIU, nanometer
from PyMieSim.utils import get_pymiescatt_sphere_dataframe
from MPSPlots.styles import mps


# Define parameters
wavelength = 632.8 * nanometer  # Wavelength of the source in meters
index = (1.4 + 0.2j) * RIU  # Refractive index of the sphere
medium_index = 1.2 * RIU  # Refractive index of the medium
optical_power = 1 * watt  # Power of the light source in watts
NA = 0.2 * AU  # Numerical aperture
diameters = np.geomspace(10, 6_000, 800) * nanometer  # Diameters from 10 nm to 6 μm

# Configure the Gaussian source
source = Gaussian(
    wavelength=wavelength,
    polarization=0 * degree,
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
experiment = Setup(scatterer=scatterer, source=source)

# Compute PyMieSim scattering efficiency data
pymiesim_dataframe = experiment.get('Qsca', 'Qext', 'Qabs', 'g', 'Qpr').reset_index('diameter')

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
        xlabel=r'Diameter [$\mu$m]',
        ylabel='Scattering Efficiency',
        title='Scattering Efficiency Comparison for Sphere Particles'
    )
    plt.legend()
    plt.show()
