"""
SPF Computation
===============

This example demonstrates the computation and visualization of the Scattering Phase Function (SPF) using PyMieSim.
"""

from PyMieSim.units import ureg

from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import RightCircular
from PyMieSim.single.representations import SPF
from PyMieSim.single import Setup

polarization = RightCircular()

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3 * ureg.AU,
)

scatterer = CoreShell(
    core_diameter=500 * ureg.nanometer,
    shell_thickness=400 * ureg.nanometer,
    core_material=1.4 * ureg.RIU,
    shell_material=1.8 * ureg.RIU,
    medium=1.0 * ureg.RIU,
)

setup = Setup(
    scatterer=scatterer,
    source=source,
)

spf = SPF(setup=setup, sampling=100)

figure = spf.plot()
