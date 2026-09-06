#!/usr/bin/env python3
"""Print the next semantic-version tag from the highest existing release tag."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
TAG_PATTERN = re.compile(r"v(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)$")


def latest_version() -> tuple[int, int, int]:
    """Return the highest semantic-version tag, including orphaned release tags."""
    completed = subprocess.run(
        ["git", "tag", "--list", "v*"],
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise RuntimeError("could not list Git tags")
    versions = []
    for tag in completed.stdout.splitlines():
        match = TAG_PATTERN.fullmatch(tag)
        if match is not None:
            versions.append(tuple(int(match.group(name)) for name in ("major", "minor", "patch")))
    if not versions:
        raise RuntimeError("no semantic-version tag was found")
    return max(versions)


def main() -> int:
    """Print the requested semantic-version increment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("major", "minor", "patch"))
    arguments = parser.parse_args()
    try:
        major, minor, patch = latest_version()
    except RuntimeError as error:
        print(f"could not determine next release: {error}", file=sys.stderr)
        return 1

    if arguments.kind == "major":
        major, minor, patch = major + 1, 0, 0
    elif arguments.kind == "minor":
        minor, patch = minor + 1, 0
    else:
        patch += 1
    print(f"v{major}.{minor}.{patch}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
