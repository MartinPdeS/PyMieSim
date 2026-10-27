"""Sphinx configuration for PyMieSim documentation."""

import subprocess
import sys
from pathlib import Path

from PyOptik.directories import user_data_path
from martinpdes_docs import configure


# Gallery examples use named PyOptik materials. Initialize their local catalog
# on fresh documentation runners before Sphinx executes those examples.
if not (user_data_path / "rii" / "catalog-nk.yml").is_file():
    subprocess.run([sys.executable, "-m", "PyOptik", "setup", "--no-progress"], check=True)


package_name = "PyMieSim"
docs_source = Path(__file__).resolve().parent
gallery_sections = [
    "single",
    "experiment",
    "inverse",
    "interpretation",
    "benchmarks",
    "validation",
    "extras",
]

globals().update(
    configure(
        project=package_name,
        package=package_name,
        version="auto",
        author="Martin Poinsinet de Sivry-Houle",
        repository=f"MartinPdeS/{package_name}",
        logo="_static/logo.png",
        favicon="_static/favicon.png",
        prefer_installed_package=True,
        gallery_sections=gallery_sections,
        gallery_default_thumb=str(docs_source / "_static" / "gallery-default.png"),
        gallery_options={"capture_repr": ()},
        logo_text="",
        logo_alt_text=f"{package_name} light-scattering logo",
        html_css_files=["layout.css"],
    )
)
