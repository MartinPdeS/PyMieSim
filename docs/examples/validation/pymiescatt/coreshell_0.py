"""
Core-Shell Particles: 0
=======================

"""

# Standard library imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PyMieSim.units import ureg

# PyMieSim imports
from PyMieSim.experiment.scatterer import CoreShellSet
from PyMieSim.experiment.source import GaussianSet, PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.directories import validation_data_path

# Define parameters
wavelength = 600 * ureg.nanometer  # Light source wavelength in meters
polarization = PolarizationSet(angles=0 * ureg.degree)
optical_power = 1 * ureg.watt  # Power in watts
NA = 0.2 * ureg.AU  # Numerical aperture
medium_index = 1.0 * ureg.RIU
core_index = 1.5 * ureg.RIU
shell_index = 1.4 * ureg.RIU
shell_thickness = 600 * ureg.nanometer  # Shell width in meters
core_diameters = np.geomspace(10, 500, 40) * ureg.nanometer  # Core diameters in meters

# Setup source
source = GaussianSet(
    wavelength=[600] * ureg.nanometer,
    polarization=polarization,
    optical_power=[1] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU
)

# Setup scatterer
scatterer = CoreShellSet(
    core_diameter=core_diameters,
    shell_thickness=[600] * ureg.nanometer,
    core_refractive_index=[1.5] * ureg.RIU,
    shell_refractive_index=[1.4] * ureg.RIU,
    medium_refractive_index=[1.0] * ureg.RIU,
    source=source,
)

# Define experiment setup
experiment = Setup(scatterer=scatterer, source=source)

comparison_measures = ["Qsca", "Qext", "Qabs", "g", "Qpr", "Qback"]

# Simulate using PyMieSim
pymiesim = experiment.get(*comparison_measures, as_numpy=True)

pymiescatt_dataframe = pd.read_csv(
    validation_data_path / "pymiescatt/example_coreshell_0.csv"
)

figure, ax = plt.subplots(1, 1)


for string in comparison_measures:
    ax.plot(
        pymiescatt_dataframe["core_diameter"],
        pymiescatt_dataframe[string],
        label="PyMieScatt: " + string,
        linewidth=3,
    )


for data, string in zip(pymiesim, comparison_measures):
    ax.plot(
        core_diameters.to(ureg.nanometer).magnitude,
        data,
        label="PyMieSim: " + string,
        linestyle="--",
        color="black",
        linewidth=1.5,
    )

ax.set(
    xlabel=r"Core Diameter [$\mu$m]",
    ylabel="Scattering Efficiency",
    title="Scattering Efficiency Comparison for Core-Shell Particles",
)
plt.legend()
plt.show()
