"""Validate the internal electric field at the origin of spherical scatterers."""

if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__]))

import numpy as np
import pytest
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.scatterer import CoreShell, Sphere
from PyMieSim.single.setup import Setup
from PyMieSim.single.source import PlaneWave
from PyMieSim.units import ureg


@pytest.mark.parametrize("scatterer_type", [Sphere, CoreShell], ids=["sphere", "core-shell"])
def test_internal_field_is_continuous_at_origin(scatterer_type):
    """The complex field and its magnitude have the same limit from several directions."""
    if scatterer_type is Sphere:
        scatterer = Sphere(diameter=400 * ureg.nanometer, material=1.6 + 0.05j, medium=1.33)
    else:
        scatterer = CoreShell(
            core_diameter=230 * ureg.nanometer, shell_thickness=85 * ureg.nanometer,
            core_material=1.7 + 0.02j, shell_material=1.6, medium=1.33,
        )
    source = PlaneWave(
        wavelength=600 * ureg.nanometer,
        polarization=PolarizationState(jones_vector=[0.6, 0.8j]),
        amplitude=2 * ureg.volt / ureg.meter,
    )
    setup = Setup(scatterer=scatterer, source=source)
    positions = np.array([
        [0., 0., 0.],
        [1e-4, 0., 0.],
        [0., -1e-4, 0.],
        [0., 0., 1e-4],
        [1e-7, -1e-7, 1e-7],
        [1e-18, 0., 0.],
    ]) * ureg.nanometer
    coordinates = dict(x=positions[:, 0], y=positions[:, 1], z=positions[:, 2])
    components = np.array([
        setup.get_total_nearfields(**coordinates, field_type=component).to("volt/meter").magnitude
        for component in ("Ex", "Ey", "Ez")
    ])
    assert np.abs(components[0, 1]) > 0.1
    assert np.abs(components[1, 1]) > 0.1
    assert components[2, 0] == 0
    for nearby in components[:, 1:].T:
        np.testing.assert_allclose(nearby, components[:, 0], rtol=0, atol=5e-6)

    magnitude = setup.get_total_nearfields(**coordinates, field_type="|E|").to("volt/meter").magnitude
    np.testing.assert_allclose(magnitude, np.linalg.norm(components, axis=0), rtol=1e-12, atol=0)
    np.testing.assert_allclose(magnitude[0], magnitude[1:], rtol=0, atol=5e-6)


def test_small_sphere_center_field_approaches_electrostatic_limit():
    """The small-sphere internal field obeys the dielectric boundary solution."""
    host_index = 1.33
    particle_index = 1.6 + 0.05j
    source = PlaneWave(
        wavelength=600 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        amplitude=2 * ureg.volt / ureg.meter,
    )
    setup = Setup(
        scatterer=Sphere(diameter=0.5 * ureg.nanometer, material=particle_index, medium=host_index),
        source=source,
    )
    origin = [0] * ureg.nanometer
    center_field = setup.get_total_nearfields(x=origin, y=origin, z=origin, field_type="Ex").to("volt/meter").magnitude[0]
    expected_field = 2 * 3 * host_index**2 / (particle_index**2 + 2 * host_index**2)
    # The leading retardation correction scales with the square of the size parameter.
    np.testing.assert_allclose(center_field, expected_field, rtol=2e-5, atol=0)
