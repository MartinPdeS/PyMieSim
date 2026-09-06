"""
Photodiode Detector
===================

This example demonstrates the initialization and visualization of a Photodiode detector using PyMieSim.
"""

from PyMieSim import (
    ureg,
    Sphere,
    Gaussian,
    PolarizationState,
    Photodiode,
    Simulation,
)



polarization_state = PolarizationState(
    angle=0 * ureg.degree,
)

source = Gaussian(
    wavelength=1550 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = Sphere(
    diameter=1800 * ureg.nanometer,
    medium=1.0,
    material=1.5,
)

detector = Photodiode(
    numerical_aperture=0.2,
    gamma_offset=0 * ureg.degree,
    phi_offset=30 * ureg.degree,
    polarization_filter=0 * ureg.degree,
    medium=3,
)

setup = Simulation(scatterer=scatterer, source=source, detector=detector)

setup.plot_system()

