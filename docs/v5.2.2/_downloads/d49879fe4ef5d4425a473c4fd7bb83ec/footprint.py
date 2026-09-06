"""
Scatterer Footprint Calculation and Visualization
=================================================

This example demonstrates how to compute and visualize the footprint of a scatterer using PyMieSim.
"""
from PyMieSim import (
    ureg,
    Sphere,
    CoherentMode,
    Gaussian,
    PolarizationState,
    Simulation,
)


polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=1 * ureg.micrometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = Sphere(
    diameter=2 * ureg.micrometer,
    medium=1.0,
    material=1.8,
)

detector = CoherentMode(
    mode_number="HG02",
    numerical_aperture=0.3,
    sampling=200,
    gamma_offset=0 * ureg.degree,
    phi_offset=0 * ureg.degree,
    rotation=0 * ureg.degree,
    medium=1.0,
)

setup = Simulation(
    scatterer=scatterer,
    source=source,
    detector=detector,
)

footprint = setup.get_representation("footprint", sampling=100)

figure = footprint.plot()
