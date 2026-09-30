"""Preflight checks for PyMieSim's compiled extension modules."""

from importlib.machinery import EXTENSION_SUFFIXES
from importlib.util import find_spec


REQUIRED_NATIVE_MODULES = (
    "PyMieSim._pint",
    "PyMieSim.coordinates",
    "PyMieSim.material",
    "PyMieSim.polarization",
    "PyMieSim.mesh",
    "PyMieSim.labeled_array",
    "PyMieSim.distributions",
    "PyMieSim.inverse",
    "PyMieSim.single.mode_field",
    "PyMieSim.single.setup",
    "PyMieSim.single.source",
    "PyMieSim.single.scatterer",
    "PyMieSim.single.optical_interface",
    "PyMieSim.single.detector",
    "PyMieSim.experiment.polarization_set",
    "PyMieSim.experiment.source_set",
    "PyMieSim.experiment.material_set",
    "PyMieSim.experiment.scatterer_set",
    "PyMieSim.experiment.detector_set",
    "PyMieSim.experiment._setup",
)

STARTUP_NATIVE_MODULES = REQUIRED_NATIVE_MODULES[:8]


def _has_extension(module_name: str) -> bool:
    """Return whether an importable compiled module is available."""
    try:
        spec = find_spec(module_name)
    except (ImportError, ModuleNotFoundError, FileNotFoundError):
        return False
    return spec is not None and any(spec.origin.endswith(suffix) for suffix in EXTENSION_SUFFIXES if spec.origin)


def missing_native_extensions(
    modules: tuple[str, ...] = REQUIRED_NATIVE_MODULES,
) -> tuple[str, ...]:
    """Return compiled modules that are not discoverable in this install."""
    return tuple(module for module in modules if not _has_extension(module))


def check_native_extensions() -> None:
    """Raise an actionable error when the native package is incomplete."""
    # Nested extension modules are imported by ``single`` and ``experiment``
    # during package initialization. Checking them here would inspect them
    # before their parent packages have finished loading.
    missing = missing_native_extensions(STARTUP_NATIVE_MODULES)
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
    "STARTUP_NATIVE_MODULES",
    "check_native_extensions",
    "missing_native_extensions",
]
