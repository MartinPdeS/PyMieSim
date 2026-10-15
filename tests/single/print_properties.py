"""The property tables use a sensible default across all scatterer geometries."""

import ctypes
import re
import sys

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
def test_print_properties_defaults_to_four_decimal_places(scatterer_type, properties, capfd):
    scatterer = scatterer_type(**properties, medium=1.0)
    scatterer.init(source=PlaneWave(
        wavelength=532 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    ))

    # Native tables use printf; flush C stdio before reading captured file descriptors.
    runtime_names = ("ucrtbase", "msvcrt") if sys.platform == "win32" else (None,)
    runtimes = [ctypes.CDLL(name) for name in runtime_names]
    for runtime in runtimes:
        runtime.fflush.argtypes = [ctypes.c_void_p]
        runtime.fflush.restype = ctypes.c_int

    def flush_native_output():
        for runtime in runtimes:
            runtime.fflush(None)

    scatterer.print_properties()
    flush_native_output()
    default_output = capfd.readouterr().out
    assert "Qsca" in default_output
    assert re.search(r"\d\.\d{4}e[+-]\d+", default_output)

    scatterer.print_properties(4)
    flush_native_output()
    assert capfd.readouterr().out == default_output

    scatterer.print_properties(precision=2)
    flush_native_output()
    shorter_output = capfd.readouterr().out
    assert shorter_output != default_output
    assert re.search(r"\d\.\d{2}e[+-]\d+", shorter_output)
