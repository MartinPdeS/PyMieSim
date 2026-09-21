"""
Integrating sphere
==================

This example demonstrates the initialization and visualization of an Integrating Sphere detector using PyMieSim.
"""
from PyMieSim import (
    ureg,
    Sphere,
    Gaussian,
    PolarizationState,
    IntegratingSphere,
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

detector = IntegratingSphere(
    sampling=200,
)

setup = Simulation(scatterer=scatterer, source=source, detector=detector)

setup.plot_system()
