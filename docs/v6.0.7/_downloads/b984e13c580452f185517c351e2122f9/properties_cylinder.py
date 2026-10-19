"""
Print cylinder properties
=========================

Initialize the scatterer with a source, then call ``print_properties()``.
The printed table uses four decimal places in scientific notation by default.
Efficiencies are dimensionless; geometric lengths and cross sections include
their units in the table. This example produces text rather than a plot.
"""

# sphinx_gallery_thumbnail_path = "_static/gallery-default.png"

from PyMieSim import (
    ureg,
    InfiniteCylinder,
    Gaussian,
    PolarizationState,
)

polarization_state = PolarizationState(angle=30 * ureg.degree)

source = Gaussian(
    wavelength=750 * ureg.nanometer,
    polarization=polarization_state,
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = InfiniteCylinder(
    diameter=300 * ureg.nanometer,
    material=(1.4 + 0.1j),
    medium=1.33,
)

scatterer.init(source=source)

scatterer.print_properties()
