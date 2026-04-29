"""
SPF Computation
===============

This example demonstrates the computation and visualization of the Scattering Phase Function (SPF) using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian

from PyMieSim.single.setup import Setup

from PyMieSim.polarization import PolarizationState
polarization = PolarizationState(angle=90 * ureg.degree)

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)


scatterer = CoreShell(
    core_diameter=300 * ureg.nanometer,
    shell_thickness=10 * ureg.nanometer,
    core_material=1.48,
    shell_material=1.38,
    medium=1.33,
)

setup = Setup(
    scatterer=scatterer,
    source=source,
)

spf = setup.get_representation("spf", sampling=100)

figure = spf.plot()
