"""
SPF Computation
===============

This example demonstrates the computation and visualization of the Scattering Phase Function (SPF) using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = CoreShell(
    core_diameter=500 * ureg.nanometer,
    shell_thickness=400 * ureg.nanometer,
    source=source,
    core_property=1.4 * ureg.RIU,
    shell_property=1.8 * ureg.RIU,
    medium_property=1.0 * ureg.RIU,
)

data = scatterer.get_spf(sampling=300)

figure = data.plot()
