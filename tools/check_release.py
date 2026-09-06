#!/usr/bin/env python3
"""Check PyMieSim's tag-derived release metadata."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys
import tomllib


ROOT = Path(__file__).resolve().parents[1]
VERSION_FILE = ROOT / "PyMieSim" / "_version.py"
TAG_PATTERN = re.compile(r"v(?P<version>(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*))$")


def normalized_version(version: str) -> str:
    """Return a version without an optional Git tag prefix."""
    return version.removeprefix("v")


def source_version() -> str:
    """Return the generated source version."""
    match = re.search(
        r"^__version__\s*=\s*version\s*=\s*['\"]([^'\"]+)",
        VERSION_FILE.read_text(encoding="utf-8"),
        flags=re.MULTILINE,
    )
    if match is None:
        raise RuntimeError("could not find the generated package version")
    return match.group(1)


def latest_tag() -> str | None:
    """Return the newest reachable semantic-version tag, if one exists."""
    completed = subprocess.run(
        ["git", "describe", "--tags", "--abbrev=0", "--match", "v[0-9]*"],
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", help="expected version, optionally prefixed with v")
    return parser.parse_args()


def main() -> int:
    """Print release metadata and return nonzero when it is inconsistent."""
    arguments = parse_arguments()
    try:
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        if "version" not in project["project"].get("dynamic", []):
            raise RuntimeError("pyproject.toml must retain dynamic SCM versioning")
        version = source_version()
    except (KeyError, RuntimeError) as error:
        print(f"release check failed: {error}", file=sys.stderr)
        return 1

    expected = normalized_version(arguments.version) if arguments.version else version
    tag = latest_tag()
    failures = []
    print(f"{'OK' if version == expected else 'MISMATCH':8} PyMieSim/_version.py: {version}")
    if version != expected:
        failures.append(f"generated version is {version}; expected {expected}")
    if tag is not None:
        tag_version = normalized_version(tag)
        print(f"{'OK' if tag_version == expected else 'MISMATCH':8} Git tag: {tag}")
        if tag_version != expected:
            failures.append(f"latest Git tag is {tag}; expected v{expected}")
    else:
        print(f"INFO     Git tag: v{expected} has not been created")

    if failures:
        print("\n".join(f"ERROR: {failure}" for failure in failures), file=sys.stderr)
        return 1
    print(f"Release metadata consistently declares {expected}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
