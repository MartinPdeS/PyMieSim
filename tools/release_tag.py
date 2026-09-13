#!/usr/bin/env python3
"""Create a PyMieSim release commit and annotated Git tag.

The tag is the canonical package version.  This helper aligns the project
metadata, Zenodo metadata, Conda recipe, and generated source version file
before committing and tagging the release.  It intentionally does not push
commits or tags to a remote.
"""

from __future__ import annotations

import argparse
from datetime import date
import json
import os
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
ZENODO_PATH = ROOT / ".zenodo.json"
CONDA_RECIPE_PATH = ROOT / "conda.recipe" / "meta.yaml"
PYPROJECT_PATH = ROOT / "pyproject.toml"
VERSION_FILE = ROOT / "PyMieSim" / "_version.py"
TAG_PATTERN = re.compile(
    r"v(?P<version>(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*))$"
)


def run(*command: str, capture_output: bool = False, env: dict[str, str] | None = None) -> str:
    """Run a repository command and return stripped standard output."""
    completed = subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=capture_output,
        env=env,
    )
    return completed.stdout.strip() if capture_output else ""


def validate_tag(tag: str) -> str:
    """Return the PEP 440 version represented by a release tag."""
    match = TAG_PATTERN.fullmatch(tag)
    if match is None:
        raise ValueError("tag must use the form vMAJOR.MINOR.PATCH, for example v5.2.2")
    return match.group("version")


def require_clean_worktree() -> None:
    """Refuse to mix a release commit with unrelated working-tree changes."""
    status = run("git", "status", "--porcelain", capture_output=True)
    if status:
        raise RuntimeError("working tree is not clean; commit or stash changes before creating a release tag")


def require_unused_tag(tag: str) -> None:
    """Refuse to overwrite an existing local tag."""
    existing = run("git", "tag", "--list", tag, capture_output=True)
    if existing:
        raise RuntimeError(f"tag {tag} already exists")


def update_zenodo(version: str) -> None:
    """Update release-specific Zenodo metadata without changing its scope."""
    metadata = json.loads(ZENODO_PATH.read_text(encoding="utf-8"))
    metadata["version"] = version
    metadata["publication_date"] = date.today().isoformat()
    ZENODO_PATH.write_text(f"{json.dumps(metadata, indent=2)}\n", encoding="utf-8")


def update_conda_recipe(version: str) -> None:
    """Write the release version into the Conda recipe."""
    recipe = CONDA_RECIPE_PATH.read_text(encoding="utf-8")
    updated_recipe, replacements = re.subn(
        r"(?m)^(?P<prefix>\s*version:\s*)\S+(?P<suffix>\s*(?:#.*)?)$",
        rf"\g<prefix>{version}\g<suffix>",
        recipe,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("could not find exactly one package version in conda.recipe/meta.yaml")
    CONDA_RECIPE_PATH.write_text(updated_recipe, encoding="utf-8")


def update_pyproject(version: str) -> None:
    """Write the release version into the Python project metadata."""
    project = PYPROJECT_PATH.read_text(encoding="utf-8")
    updated_project, replacements = re.subn(
        r'(?ms)(^\[project\]\n.*?^version\s*=\s*)["\'][^"\']+["\']',
        rf'\g<1>"{version}"',
        project,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("could not find exactly one project version in pyproject.toml")
    PYPROJECT_PATH.write_text(updated_project, encoding="utf-8")


def generate_version_file(version: str) -> None:
    """Generate ``_version.py`` through the configured SCM-versioning tool."""
    environment = os.environ.copy()
    environment["SETUPTOOLS_SCM_PRETEND_VERSION"] = version
    run(
        sys.executable,
        "-m",
        "vcs_versioning",
        "--force-write-version-files",
        env=environment,
    )
    if not VERSION_FILE.exists():
        raise RuntimeError("SCM versioning did not generate PyMieSim/_version.py")


def create_release(tag: str, version: str) -> None:
    """Write metadata, commit it, and create the annotated release tag."""
    update_zenodo(version)
    update_conda_recipe(version)
    update_pyproject(version)
    generate_version_file(version)
    run(
        "git",
        "add",
        ".zenodo.json",
        "conda.recipe/meta.yaml",
        "pyproject.toml",
        "PyMieSim/_version.py",
    )
    run("git", "commit", "-m", f"Release {tag}")
    run("git", "tag", "-a", tag, "-m", f"Release {tag}")


def parse_arguments() -> argparse.Namespace:
    """Parse the requested release tag."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="annotated release tag, for example v5.2.2")
    return parser.parse_args()


def main() -> int:
    """Create the requested release tag after conservative preflight checks."""
    arguments = parse_arguments()
    try:
        version = validate_tag(arguments.tag)
        require_clean_worktree()
        require_unused_tag(arguments.tag)
        create_release(arguments.tag, version)
    except (RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"release aborted: {error}", file=sys.stderr)
        return 1

    print(f"created release commit and annotated tag {arguments.tag}")
    print("Push it when ready with: git push origin HEAD --tags")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
