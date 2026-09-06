#!/usr/bin/env python3
"""Create a PyMieSim release commit and annotated semantic-version tag."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
VERSION_FILE = ROOT / "PyMieSim" / "_version.py"
TAG_PATTERN = re.compile(r"v(?P<version>(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*))$")


def run(*command: str, capture_output: bool = False, env: dict[str, str] | None = None) -> str:
    """Run a repository command and return stripped standard output."""
    completed = subprocess.run(
        command, cwd=ROOT, check=True, text=True, capture_output=capture_output, env=env
    )
    return completed.stdout.strip() if capture_output else ""


def validate_tag(tag: str) -> str:
    """Return the PEP 440 version represented by a release tag."""
    match = TAG_PATTERN.fullmatch(tag)
    if match is None:
        raise ValueError("tag must use vMAJOR.MINOR.PATCH, for example v5.2.1")
    return match.group("version")


def require_clean_worktree() -> None:
    """Refuse to mix a release commit with unrelated changes."""
    if run("git", "status", "--porcelain", capture_output=True):
        raise RuntimeError("working tree is not clean; commit or stash changes before creating a release tag")


def require_unused_tag(tag: str) -> None:
    """Refuse to overwrite an existing local tag."""
    if run("git", "tag", "--list", tag, capture_output=True):
        raise RuntimeError(f"tag {tag} already exists")


def generate_version_file(version: str) -> None:
    """Generate the tracked source version through SCM versioning."""
    environment = os.environ.copy()
    environment["SETUPTOOLS_SCM_PRETEND_VERSION"] = version
    run(sys.executable, "-m", "vcs_versioning", "--force-write-version-files", env=environment)
    if not VERSION_FILE.exists():
        raise RuntimeError("SCM versioning did not generate PyMieSim/_version.py")


def main() -> int:
    """Create a release commit and annotated tag without pushing either."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="annotated release tag, for example v5.2.1")
    arguments = parser.parse_args()
    try:
        version = validate_tag(arguments.tag)
        require_clean_worktree()
        require_unused_tag(arguments.tag)
        generate_version_file(version)
        run("git", "add", "PyMieSim/_version.py")
        run("git", "commit", "-m", f"Release {arguments.tag}")
        run("git", "tag", "-a", arguments.tag, "-m", f"Release {arguments.tag}")
    except (RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"release aborted: {error}", file=sys.stderr)
        return 1
    print(f"created release commit and annotated tag {arguments.tag}")
    print("Push it when ready with: git push origin HEAD --tags")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
