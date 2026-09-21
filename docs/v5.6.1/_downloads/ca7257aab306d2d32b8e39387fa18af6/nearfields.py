"""
Near-Field Computation and Visualization
=========================================
"""

from PyMieSim import (
    ureg,
    PolarizationState,
    CoreShell,
    Gaussian,
    Simulation,
)


source = Gaussian(
    wavelength=100 * ureg.nanometer,
    polarization=PolarizationState(angle=0 * ureg.degree),
    optical_power=1 * ureg.watt,
    numerical_aperture=0.3,
)

scatterer = CoreShell(
    core_diameter=200 * ureg.nanometer,
    shell_thickness=100 * ureg.nanometer,
    core_material=1.3 + 2j,
    shell_material=1.5,
    medium=1.0,
)

setup = Simulation(
    scatterer=scatterer,
    source=source,
)

near_field = setup.get_representation("nearfields")

near_field.plot(
    "Ez:real",
    "Ex:real",
    type="total",
    plane_origin=(0.0, 0.0, 0.0),
    plane_normal=(0.0, 1.0, 0.0),
    sampling=100,
    extent_scale=2,
    tight_layout=True,
)