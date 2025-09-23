"""
Samples Properties
==================
"""

from TypedUnit import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = Sphere(
    diameter=800 * ureg.nanometer,
    source=source,
    property=1.4 * ureg.RIU,
    medium_property=1.0 * ureg.RIU,
)

# %%
# Plotting the farfield pattern
farfield = scatterer.get_farfield(sampling=300)
farfield.plot()

# %%
# Printing the scatterer properties
scatterer.print_properties()
