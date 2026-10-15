"""The property tables use a sensible default across all scatterer geometries."""

import re

import pytest

from PyMieSim import CoreShell, InfiniteCylinder, PlaneWave, PolarizationState, Sphere, ureg


@pytest.mark.parametrize(
    ("scatterer_type", "properties"),
    [
        (Sphere, {"diameter": 200 * ureg.nanometer, "material": 1.5}),
        (InfiniteCylinder, {"diameter": 200 * ureg.nanometer, "material": 1.5}),
        (CoreShell, {
            "core_diameter": 100 * ureg.nanometer,
            "shell_thickness": 50 * ureg.nanometer,
            "core_material": 1.5,
            "shell_material": 1.4,
        }),
    ],
)
def test_print_properties_defaults_to_four_decimal_places(scatterer_type, properties, capsys):
    scatterer = scatterer_type(**properties, medium=1.0)
    scatterer.init(source=PlaneWave(
        wavelength=532 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    ))

    scatterer.print_properties()
    default_output = capsys.readouterr().out
    assert "Qsca" in default_output
    assert re.search(r"\d\.\d{4}e[+-]\d+", default_output)

    scatterer.print_properties(4)
    assert capsys.readouterr().out == default_output

    scatterer.print_properties(precision=2)
    shorter_output = capsys.readouterr().out
    assert shorter_output != default_output
    assert re.search(r"\d\.\d{2}e[+-]\d+", shorter_output)
