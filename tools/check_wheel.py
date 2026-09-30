#!/usr/bin/env python3
"""Verify that a PyMieSim wheel contains every native extension."""

import argparse
from pathlib import Path
import zipfile


NATIVE_MODULES = (
    "PyMieSim/_pint",
    "PyMieSim/coordinates",
    "PyMieSim/material",
    "PyMieSim/polarization",
    "PyMieSim/mesh",
    "PyMieSim/labeled_array",
    "PyMieSim/distributions",
    "PyMieSim/inverse",
    "PyMieSim/single/mode_field",
    "PyMieSim/single/setup",
    "PyMieSim/single/source",
    "PyMieSim/single/scatterer",
    "PyMieSim/single/optical_interface",
    "PyMieSim/single/detector",
    "PyMieSim/experiment/polarization_set",
    "PyMieSim/experiment/source_set",
    "PyMieSim/experiment/material_set",
    "PyMieSim/experiment/scatterer_set",
    "PyMieSim/experiment/detector_set",
    "PyMieSim/experiment/_setup",
)


def missing_modules(wheel: Path) -> list[str]:
    """Return native module names absent from a wheel archive."""
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
    return [
        module
        for module in NATIVE_MODULES
        if not any(
            name.startswith(f"{module}.")
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
