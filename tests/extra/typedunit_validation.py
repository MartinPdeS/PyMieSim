"""TypedUnit validation precedes numerical conversion in native entry points."""

import numpy as np
import pint
import pytest
import TypedUnit
from PyMieSim import ConstantMaterial, ConstantMedium, PlaneWave, PolarizationState, Simulation, Sphere, ureg
from PyMieSim.single.source import Gaussian


def test_registry_and_equivalent_units_preserve_scattering():
    """Equivalent length units preserve both inputs and the computed response."""
    assert ureg is TypedUnit.ureg
    diameter = 0.1 * ureg.micrometer
    wavelength = 532 * ureg.nanometer
    source = PlaneWave(
        wavelength=wavelength,
        polarization=PolarizationState(angle=0 * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    )
    first_sphere = Sphere(diameter=diameter, material=1.5 + 0.01j, medium=1.0)
    second_sphere = Sphere(diameter=100 * ureg.nanometer, material=1.5 + 0.01j, medium=1.0)
    first_simulation = Simulation(scatterer=first_sphere, source=source)
    second_simulation = Simulation(scatterer=second_sphere, source=source)
    first_result = first_simulation.run("Qsca")
    second_result = second_simulation.run("Qsca")
    np.testing.assert_allclose(first_result.magnitude, second_result.magnitude)
    assert diameter.units == ureg.micrometer
    assert wavelength.units == ureg.nanometer


@pytest.mark.parametrize(
    "diameter",
    [100, np.nan * ureg.nanometer, np.array([100]) * ureg.nanometer, 1 * ureg.second, 100 * pint.UnitRegistry().nanometer],
)
def test_sphere_validates_dimensional_scalar_before_conversion(diameter):
    """Lengths require a finite scalar quantity from the shared registry."""
    with pytest.raises(ValueError):
        Sphere(diameter=diameter, material=1.5, medium=1.0)


@pytest.mark.parametrize("constructor", [ConstantMaterial, ConstantMedium])
@pytest.mark.parametrize("value", [True, np.nan, [1.5], 1.5 * ureg.dimensionless, 1.5 * pint.UnitRegistry().dimensionless])
def test_constant_indices_require_plain_finite_scalars(constructor, value):
    """All Pint registries are rejected for deliberately unitless indices."""
    with pytest.raises(ValueError, match="refractive_index"):
        constructor(refractive_index=value)


@pytest.mark.parametrize("angle", [0, 0 * ureg.dimensionless, np.nan * ureg.degree])
def test_polarization_requires_explicit_finite_angles(angle):
    """Pint dimensionlessness cannot stand in for explicit angular units."""
    with pytest.raises(ValueError):
        PolarizationState(angle=angle)


@pytest.mark.parametrize("numerical_aperture", [True, np.nan, 0.2 * ureg.dimensionless])
def test_gaussian_aperture_is_a_plain_finite_scalar(numerical_aperture):
    """Apertures are checked before casting to native doubles."""
    with pytest.raises(ValueError, match="numerical_aperture"):
        Gaussian(
            wavelength=532 * ureg.nanometer,
            polarization=PolarizationState(angle=0 * ureg.degree),
            optical_power=1 * ureg.watt,
            numerical_aperture=numerical_aperture,
        )
