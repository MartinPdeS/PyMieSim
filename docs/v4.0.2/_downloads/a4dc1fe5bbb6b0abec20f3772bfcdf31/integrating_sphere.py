"""
Integrating sphere
==================

This example demonstrates the initialization and visualization of an Integrating Sphere detector using PyMieSim.
"""
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single import Setup


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

setup = Setup(scatterer=scatterer, source=source, detector=detector)

setup.plot_system()
