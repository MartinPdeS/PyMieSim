"""
Near-Field Computation and Visualization
=========================================
"""

from PyMieSim.units import ureg
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.setup import Setup


source = Gaussian(
    wavelength=100 * ureg.nanometer,
    polarization=PolarizationState(angle=0 * ureg.degree),
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = InfiniteCylinder(
    diameter=200 * ureg.nanometer,
    material=1.8 + 2j,
    medium=1.0,
)

setup = Setup(
    scatterer=scatterer,
    source=source,
)

near_field = setup.get_representation("nearfields")

near_field.plot(
    "|E|:abs",
    type="total",
    plane_origin=(0.0, 0.0, 0.0),
    plane_normal=(1.0, 1.0, 0.0),
    sampling=100,
    extent_scale=2,
    tight_layout=True,
)