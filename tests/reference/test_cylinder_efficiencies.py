"""Independent partial-wave reference for a cylinder at normal incidence."""

import numpy as np
import pytest
from scipy.special import hankel1, h1vp, jv, jvp
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.setup import Setup
from PyMieSim.single.source import PlaneWave
from PyMieSim.units import ureg


@pytest.mark.parametrize("angle", [0, 90], ids=["transverse-E", "axial-E"])
def test_cylinder_scattering_and_extinction(angle):
    """Match the dimensionless scattering and extinction widths divided by diameter."""
    diameter, wavelength = 300., 600.
    host, particle = 1.33, 1.6 + 0.05j
    setup = Setup(
        scatterer=InfiniteCylinder(diameter=diameter * ureg.nanometer, material=particle, medium=host),
        source=PlaneWave(wavelength=wavelength * ureg.nanometer,
                         polarization=PolarizationState(angle=angle * ureg.degree),
                         amplitude=1 * ureg.volt / ureg.meter),
    )
    x = np.pi * diameter * host / wavelength
    m = particle / host
    # Signed orders include the +/-n degeneracy. At x~2 the tails are negligible.
    orders = np.arange(-20, 21)
    interior_derivative = jvp(orders, m*x) / jv(orders, m*x)
    boundary_derivative = interior_derivative * (1/m if angle == 0 else m)
    coefficients = (jvp(orders, x) - boundary_derivative * jv(orders, x)) / (
        h1vp(orders, x) - boundary_derivative * hankel1(orders, x)
    )
    np.testing.assert_allclose(setup.get("Qsca"), 2/x * np.sum(abs(coefficients)**2), rtol=1e-10, atol=0)
    np.testing.assert_allclose(setup.get("Qext"), 2/x * coefficients.real.sum(), rtol=1e-10, atol=0)
