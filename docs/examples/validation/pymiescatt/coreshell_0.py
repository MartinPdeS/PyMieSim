"""
PyMieSim vs PyMieScatt for Core-Shell Particles
===============================================
"""

# Standard library imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import degree, watt, AU, RIU, nanometer
from PyMieSim.directories import validation_data_path
from MPSPlots.styles import mps

# Define parameters
wavelength = 600 * nanometer  # Light source wavelength in meters
polarization = 0 * degree
optical_power = 1 * watt  # Power in watts
NA = 0.2 * AU  # Numerical aperture
medium_index = 1.0 * RIU
core_index = 1.5 * RIU
shell_index = 1.4 * RIU
shell_width = 600 * nanometer  # Shell width in meters
core_diameters = np.geomspace(10, 500, 40) * nanometer  # Core diameters in meters

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
    core_property=core_index,
    shell_property=shell_index,
    medium_property=medium_index,
    source=source
)

# Define experiment setup
experiment = Setup(scatterer=scatterer, source=source)

comparison_measures = ['Qsca', 'Qext', 'Qabs', 'g', 'Qpr', 'Qback']

# Simulate using PyMieSim
pymiesim_dataframe = experiment.get(*comparison_measures).pint.dequantify().reset_index().pint.quantify()

pymiescatt_dataframe = pd.read_csv(validation_data_path / 'pymiescatt/example_coreshell_0.csv')

# Plot results
with plt.style.context(mps):
    figure, ax = plt.subplots(1, 1)


pymiescatt_dataframe.plot(x='core_diameter', y=comparison_measures, ax=ax, linewidth=3)
pymiesim_dataframe.plot(x='scatterer:core_diameter', ax=ax, color='black', linestyle='--', linewidth=1.5)

ax.set(
    xlabel=r'Core Diameter [$\mu$m]',
    ylabel='Scattering Efficiency',
    title='Scattering Efficiency Comparison for Core-Shell Particles'
)
plt.legend()
plt.show()
