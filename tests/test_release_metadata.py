"""Keep explicit release metadata consistent without Git version discovery."""

import json
import re
import runpy
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[1]


def test_explicit_release_metadata_agrees():
    """Source imports, distribution metadata, and Conda declare one version."""
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text())
    project = metadata["project"]
    assert "version" not in project.get("dynamic", [])
    version_file = ROOT / metadata["tool"]["mps-maintenance"]["version-file"]
    package = runpy.run_path(str(version_file))
    assert package["__version__"] == package["version"] == project["version"]
    assert package["version_tuple"] == tuple(int(part) for part in project["version"].split("."))
    recipe = (ROOT / "conda.recipe/meta.yaml").read_text()
    recipe_version = re.search(r"^  version: (.+)$", recipe, re.MULTILINE)[1].strip('"\'')
    assert recipe_version == project["version"]
    assert metadata["tool"]["mps-maintenance"]["version-method"] == "literal"
    dependencies = metadata["build-system"]["requires"] + project["dependencies"]
    assert not any("setuptools-scm" in item.replace("_", "-") for item in dependencies)

    zenodo = json.loads((ROOT / ".zenodo.json").read_text())
    assert zenodo["version"] == project["version"]
