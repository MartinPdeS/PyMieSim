"""
Print sphere properties
=======================

Initialize the scatterer with a source, then call ``print_properties()``.
The printed table uses four decimal places in scientific notation by default.
Efficiencies are dimensionless; geometric lengths and cross sections include
their units in the table. This example produces text rather than a plot.
"""

# sphinx_gallery_thumbnail_path = "_static/gallery-default.png"

from PyMieSim import (
    ureg,
    Sphere,
    Gaussian,
    PolarizationState,
)


polarization_state = PolarizationState(angle=0 * ureg.degree)

source = Gaussian(
    wavelength=1000 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = Sphere(
    diameter=800 * ureg.nanometer,
    material=1.4,
    medium=1.0,
)

scatterer.init(source=source)

scatterer.print_properties()
