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
        author="Martin Poinsinet de Sivry-Houle",
        repository="MartinPdeS/PyMieSim",
        logo="_static/logo.png",
        favicon="_static/favicon.png",
        prefer_installed_package=True,
        gallery_sections=gallery_sections,
        # Sphinx-Gallery resolves this path from the process working directory
        # (the repository root in CI), so pass an absolute source-tree path.
        gallery_default_thumb=str(docs_source / "_static" / "gallery-default.png"),
        logo_alt_text="PyMieSim light-scattering logo",
        autodoc_default_options={
            "members": True,
            "exclude-members": "__annotations__",
            "members-order": "bysource",
        },
        html_css_files=["layout.css"],
        html_theme_options={
            "show_nav_level": 1,
            "navbar_end": ["version-switcher", "theme-switcher", "navbar-icon-links"],
            "show_prev_next": True,
        },
    )
)

autodoc_mock_imports = ["numpy", "matplotlib", "numpydoc"]
autodoc_typehints = "description"
napoleon_use_ivar = True
show_authors = True

intersphinx_mapping = {
    "numpy": ("https://numpy.org/devdocs", None),
}
