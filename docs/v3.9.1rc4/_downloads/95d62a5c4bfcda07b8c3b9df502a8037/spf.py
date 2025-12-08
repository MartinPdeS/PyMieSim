"""
SPF Computation
===============

This example demonstrates the computation and visualization of the Scattering Phase Function (SPF) using PyMieSim.
"""

# %%
# Importing the package: PyMieSim
from TypedUnit import ureg

from PyMieSim.single.scatterer import CoreShell
from PyMieSim.single.source import Gaussian

source = Gaussian(
    wavelength=1000 * ureg.nanometer,  # 1000 nm
    polarization=0 * ureg.degree,  # Linear polarization angle in radians
    optical_power=1 * ureg.watt,  # Arbitrary units
    NA=0.3 * ureg.AU,  # Numerical Aperture
)

scatterer = CoreShell(
    core_diameter=500 * ureg.nanometer,  # 500 nm
    shell_thickness=100 * ureg.nanometer,  # 100 nm
    source=source,
    core_property=1.4 * ureg.RIU,  # Refractive property of the core
    shell_property=1.8 * ureg.RIU,  # Refractive property of the shell
    medium_property=1.0 * ureg.RIU,  # Refractive property of the surrounding medium
)

data = scatterer.get_spf(sampling=300)  # Specify the number of sampling points

figure = data.plot()
