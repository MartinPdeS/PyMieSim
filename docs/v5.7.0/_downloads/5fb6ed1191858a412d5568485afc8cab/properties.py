"""
Far-Fields Computation and Visualization
========================================

This example demonstrates the process of computing and visualizing the far-fields of a scatterer using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim import (
    ureg,
    Gaussian,
    PolarizationState,
    Sphere,
    Simulation,
)

polarization = PolarizationState(angle=30 * ureg.degree)

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = Sphere(
    diameter=1500 * ureg.nanometer,
    material=1.4,
    medium=1.0,
)

setup = Simulation(
    scatterer=scatterer,
    source=source,
)

for property in ["Qext", "Qsca", "Qabs", "Qback", "g", "Cext", "Csca", "Cabs", "Cback"]:
    value = setup.get(property)
    print(f"{property}: {value}")


