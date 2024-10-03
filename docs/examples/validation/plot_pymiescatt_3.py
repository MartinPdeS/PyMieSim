"""
PyMieSim vs PyMieScatt Scattering Efficiency Comparison for Core-Shell Particles
================================================================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import degree, watt, AU, RIU, nanometer
from PyMieSim.utils import get_pymiescatt_coreshell_dataframe
from MPSPlots.styles import mps

# Define parameters
wavelength = 600 * nanometer  # Light source wavelength in meters
polarization = 0 * degree
optical_power = 1 * watt  # Optical power in watts
NA = 0.2 * AU  # Numerical aperture
medium_index = 1.0 * RIU
core_index = 1.5 * RIU
shell_index = 1.4 * RIU
shell_width = 600 * nanometer  # Shell width in meters
core_diameters = np.geomspace(10, 500, 400) * nanometer  # Core diameters in meters

# Configure the Gaussian source
source = Gaussian(
    wavelength=wavelength,
    polarization=polarization,
    optical_power=optical_power,
    NA=NA
)

# Setup core-shell scatterer
scatterer = CoreShell(
    core_diameter=core_diameters,
    shell_width=shell_width,
    core_property=core_index,
    shell_property=shell_index,
    medium_property=medium_index,
    source=source
)

# Create experimental setup
experiment = Setup(scatterer=scatterer, source=source)

# Simulate using PyMieSim
pymiesim_dataframe = experiment.get('Qsca', 'Qext', 'Qabs', 'Qpr', 'g', 'Qback').reset_index('core_diameter')

pymiescatt_dataframe = get_pymiescatt_coreshell_dataframe(
    wavelengths=wavelength,
    core_indexes=core_index,
    shell_indexes=shell_index,
    core_diameters=core_diameters,
    shell_widths=shell_width,
    medium_indexes=medium_index
).reset_index('core_diameter')

# Plot results
with plt.style.context(mps):
    figure, ax = plt.subplots(1, 1)
    pymiescatt_dataframe.plot(x=0, ax=ax, linewidth=3)
    pymiesim_dataframe.plot(x=0, ax=ax, color='black', linestyle='--', linewidth=1.5)
    ax.set(
        xlabel='Core Diameter (μm)',
        ylabel='Scattering Efficiency',
        title='Comparison of Scattering Efficiency for Core-Shell Particles',
    )
    plt.legend()
    plt.show()
