"""
Sphere Particles: 0
===================

"""

# Standard library imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from TypedUnit import ureg

# PyMieSim imports
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.directories import validation_data_path
from MPSPlots.styles import mps


# Define parameters
wavelength = 632.8 * ureg.nanometer  # Wavelength of the light source in meters
polarization_value = 0 * ureg.degree
optical_power = 1e-3 * ureg.watt  # Power in watts
NA = 0.2 * ureg.AU  # Numerical aperture
medium_index = 1.21 * ureg.RIU
index = 1.4 * ureg.RIU
diameters = (
    np.geomspace(10, 1_000, 50) * ureg.nanometer
)  # Geometric space for diameters

# Setup source
source = Gaussian(
    wavelength=wavelength,
    polarization=polarization_value,
    optical_power=optical_power,
    NA=NA,
)

# Setup scatterer
scatterer = Sphere(
    diameter=diameters, property=index, medium_property=medium_index, source=source
)

# Define experiment setup
experiment = Setup(scatterer=scatterer, source=source)

comparison_measures = ["Qsca", "Qext", "Qabs", "g", "Qpr", "Qback"]

# Simulate using PyMieSim
pymiesim_dataframe = experiment.get(*comparison_measures)

pymiescatt_dataframe = pd.read_csv(
    validation_data_path / "pymiescatt/example_shpere_0.csv"
)

# Plot results
with plt.style.context(mps):
    figure, ax = plt.subplots(1, 1)

pymiescatt_dataframe.diameter *= 1e9

pymiescatt_dataframe.plot(x="diameter", y=comparison_measures, ax=ax, linewidth=3)
pymiesim_dataframe.plot(
    x="scatterer:diameter",
    axes=ax,
    color="black",
    linestyle="--",
    linewidth=1.5,
    show=False,
)

ax.set(
    xlabel="Diameter [nm]",
    ylabel="Scattering Efficiencies",
    title="Scattering parameters Comparison for Sphere Particles",
)

plt.legend()
plt.show()
