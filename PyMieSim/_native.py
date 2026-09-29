"""Preflight checks for PyMieSim's compiled extension modules."""

from importlib.util import find_spec


REQUIRED_NATIVE_MODULES = (
    "_pint",
    "coordinates",
    "material",
    "polarization",
    "mesh",
    "labeled_array",
    "distributions",
    "inverse",
    "mode_field",
    "setup_single",
    "source",
    "scatterer",
    "optical_interface",
    "detector",
    "polarization_set",
    "source_set",
    "material_set",
    "scatterer_set",
    "detector_set",
    "_setup",
)


def missing_native_extensions() -> tuple[str, ...]:
    """Return compiled modules that are not discoverable in this install."""
    return tuple(
        f"PyMieSim.{module}"
        for module in REQUIRED_NATIVE_MODULES
        if find_spec(f"PyMieSim.{module}") is None
    )


def check_native_extensions() -> None:
    """Raise an actionable error when the native package is incomplete."""
    missing = missing_native_extensions()
    if not missing:
        return

    missing_list = ", ".join(missing)
    raise ImportError(
        "PyMieSim's compiled native extensions are missing: "
        f"{missing_list}. This usually means that PyMieSim was imported "
        "from a source checkout before it was built, or that the package "
        "was built for a different Python interpreter. Install a released "
        "wheel with `python -m pip install PyMieSim`, or from this checkout "
        "run `make editable` with the same Python used to run your code. "
        "For a source build, see the 'Developer setup' section in "
        "README.rst."
    )


__all__ = [
    "REQUIRED_NATIVE_MODULES",
    "check_native_extensions",
    "missing_native_extensions",
]
