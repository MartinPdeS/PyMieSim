#!/usr/bin/env python3
"""Verify that a PyMieSim wheel contains every native extension."""

import argparse
from pathlib import Path
import zipfile


NATIVE_MODULES = (
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


def missing_modules(wheel: Path) -> list[str]:
    """Return native module names absent from a wheel archive."""
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
    return [
        module
        for module in NATIVE_MODULES
        if not any(
            name.startswith(f"PyMieSim/{module}.")
            and name.endswith((".so", ".pyd", ".dylib"))
            for name in names
        )
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    arguments = parser.parse_args()

    missing = missing_modules(arguments.wheel)
    if missing:
        print(f"{arguments.wheel}: missing native modules: {', '.join(missing)}")
        return 1
    print(f"{arguments.wheel}: all {len(NATIVE_MODULES)} native modules are present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
