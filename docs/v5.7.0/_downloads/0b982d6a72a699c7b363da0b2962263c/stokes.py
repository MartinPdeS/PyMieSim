"""
Stokes Parameters Computation
=============================

This example demonstrates the computation and visualization of the Stokes parameters using PyMieSim.
"""

from PyMieSim import (
    ureg,
    Sphere,
    Gaussian,
    PolarizationState,
    Simulation,
)


polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=750 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = Sphere(
    diameter=600 * ureg.nanometer,
    medium=1.0,
    material=1.4,
)

setup = Simulation(
    scatterer=scatterer,
    source=source,
)

stokes = setup.get_representation("stokes", sampling=100)

figure = stokes.plot()
