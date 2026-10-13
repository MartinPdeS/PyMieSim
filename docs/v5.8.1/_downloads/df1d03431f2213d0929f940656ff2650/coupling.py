"""
Coupling
========

This example demonstrates the process of computing the coupling efficiency of a scatterer to a detector using PyMieSim.
"""

from PyMieSim import (
    ureg,
    Gaussian,
    PolarizationState,
    Sphere,
    Photodiode,
    IntegratingSphere,
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
    material=1.45,
    medium=1.0,
)





# First detector
side_detector = Photodiode(
    numerical_aperture=0.3,
    phi_offset=90 * ureg.degree,
    gamma_offset=0 * ureg.degree,
    sampling=1000
)

setup = Simulation(
    scatterer=scatterer,
    source=source,
    detector=side_detector
)

side_coupling = setup.get("coupling")


# Second detector
forward_detector = IntegratingSphere(
    sampling=1000
)

setup = Simulation(
    scatterer=scatterer,
    source=source,
    detector=forward_detector
)

forward_coupling = setup.get("coupling")

print(f"Side coupling: {side_coupling} \n Forward coupling: {forward_coupling}  \n Ratio: {side_coupling / forward_coupling * 100} %")
