"""Validate directional scattering efficiencies and signed angular moments."""

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


def angular_asymmetry(*, setup):
    """Integrate the normalized first angular moment in cos(scattering angle)."""
    cosine, weights = np.polynomial.legendre.leggauss(64)
    perpendicular, parallel = setup.get_s1s2(angles=np.arccos(cosine) * ureg.radian)
    intensity = np.abs(perpendicular.magnitude)**2 + np.abs(parallel.magnitude)**2
    return np.dot(weights, intensity * cosine) / np.dot(weights, intensity)


@pytest.mark.parametrize("scatterer_type", [Sphere, CoreShell], ids=["sphere", "core-shell"])
@pytest.mark.parametrize("max_order", [1, 2])
def test_directional_efficiencies_and_asymmetry_include_last_retained_order(scatterer_type, max_order):
    """Compare the truncated coefficient sums with the corresponding angular amplitudes."""
    source = PlaneWave(wavelength=600 * ureg.nanometer, polarization=PolarizationState(angle=25 * ureg.degree),
                       amplitude=1 * ureg.volt / ureg.meter)
    if scatterer_type is Sphere:
        scatterer = Sphere(diameter=400 * ureg.nanometer, material=1.6, medium=1.33)
    else:
        scatterer = CoreShell(
            core_diameter=230 * ureg.nanometer, shell_thickness=85 * ureg.nanometer,
            core_material=1.7, shell_material=1.6, medium=1.33,
        )
    setup = Setup(scatterer=scatterer, source=source)
    scatterer.init(source=source, max_order=max_order)

    perpendicular, _ = setup.get_s1s2(angles=[0, 180] * ureg.degree)
    size_parameter = setup.get("size_parameter").magnitude
    expected_efficiencies = 4 * np.abs(perpendicular.magnitude)**2 / size_parameter**2
    actual_efficiencies = [setup.get(measure).magnitude for measure in ("Qforward", "Qback")]
    np.testing.assert_allclose(actual_efficiencies, expected_efficiencies, rtol=1e-12, atol=0)
    np.testing.assert_allclose(setup.get("g").magnitude, angular_asymmetry(setup=setup), rtol=1e-12, atol=0)


@pytest.mark.parametrize(("material", "expected_sign"), [(3.0, -1), (1.6, 1)])
def test_farfield_asymmetry_uses_signed_cosine_from_forward_axis(material, expected_sign):
    """Check the signed radiation moment for backward- and forward-biased scattering."""
    setup = Setup(
        scatterer=Sphere(diameter=350 * ureg.nanometer, material=material, medium=1.0),
        source=PlaneWave(
            wavelength=600 * ureg.nanometer, polarization=PolarizationState(angle=0 * ureg.degree),
            amplitude=1 * ureg.volt / ureg.meter,
        ),
    )

    expected = angular_asymmetry(setup=setup)
    assert np.sign(expected) == expected_sign
    # The sampled far-field integral uses a finite uniform elevation grid.
    np.testing.assert_allclose(setup.get("g_with_farfields").magnitude, expected, rtol=0, atol=5e-6)
