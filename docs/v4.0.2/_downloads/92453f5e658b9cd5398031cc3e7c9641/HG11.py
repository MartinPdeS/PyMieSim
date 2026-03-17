"""
Hermite-Gauss 31 Mode Detector
==============================

This example demonstrates the initialization and visualization of HG31 Mode detector using PyMieSim.
"""
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.detector import CoherentMode
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

detector = CoherentMode(
    mode_number="HG11",
    numerical_aperture=0.2,
    gamma_offset=0 * ureg.degree,
    phi_offset=30 * ureg.degree,
    rotation=0 * ureg.degree,
    polarization_filter=0 * ureg.degree,
    medium=1.0,
)

setup = Setup(scatterer=scatterer, source=source, detector=detector)

setup.plot_system()
