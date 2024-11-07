#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from MPSPlots.styles import use_mpsplots_style
import PyMieSim
from pathlib import Path
from PyMieSim.directories import doc_css_path

sys.path.append(str(Path(".").resolve()))


def setup(app):
    app.add_css_file(str(doc_css_path))


autodoc_mock_imports = [
    'numpy',
    'matplotlib',
    'MPSPlots',
]

project = 'PyMieSim'
copyright = '2021, Martin Poinsinet de Sivry-Houle'
author = 'Martin Poinsinet de Sivry-Houle'
today_fmt = '%B %d, %Y'

version = PyMieSim.__version__

extensions = [
    'sphinx.ext.mathjax',
    'pyvista.ext.plot_directive',
    'sphinx_gallery.gen_gallery',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
]

# Napoleon settings for docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True

html_logo = "_static/thumbnail.png"
html_favicon = "_static/thumbnail.png"


def reset_mpl(gallery_conf, fname):
    use_mpsplots_style()


try:
    import pyvista
    if sys.platform in ["linux", "linux2"]:
        pyvista.start_xvfb()  # Works only on linux system!
except ImportError:
    print('Could not load pyvista library for 3D rendering')


examples_files = [
    'single', 'experiment', 'validation', 'extras'
]

sphinx_gallery_conf = {
    "examples_dirs": ['../examples/' + f for f in examples_files],
    "gallery_dirs": ['gallery/' + f for f in examples_files],
    'image_scrapers': ('matplotlib', 'pyvista'),
    'filename_pattern': r'.*\.py',
    'ignore_pattern': '/__',
    'plot_gallery': True,
    'reset_modules': reset_mpl,
    'thumbnail_size': [600, 600],
    'download_all_examples': False,
    'line_numbers': False,
    'remove_config_comments': True,
    'capture_repr': ('_repr_html_', '__repr__'),
    'nested_sections': True,
}

autodoc_default_options = {
    'members': False,
    'members-order': 'bysource',
    'undoc-members': False,
    'show-inheritance': True,
    'ignore-module-all': True
}

autosectionlabel_prefix_document = True
numpydoc_show_class_members = False
add_module_names = False

source_suffix = '.rst'
master_doc = 'index'
language = 'en'
highlight_language = 'python3'
html_theme = "pydata_sphinx_theme"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

exclude_trees = []
# default_role = "autolink"
pygments_style = "sphinx"

# -- Sphinx-gallery configuration --------------------------------------------
binder_branch = "master"

major, minor = version[:2]
binder_branch = f"v{major}.{minor}.x"

html_theme_options = dict()

html_theme_options['logo'] = dict(text='PyMieSim', image="_static/logo-dark.svg")
html_theme_options["show_nav_level"] = 0

html_theme_options.update({
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/MartinPdeS/PyMieSim",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/pymiesim/",
            "icon": "fa-solid fa-box",
        },
        {
            "name": "Anaconda",
            "url": "https://anaconda.org/MartinPdeS/pymiesim",
            "icon": "fa-brands fa-python",
        },
    ],
    "navbar_align": "left",
    "navbar_end": ["version-switcher", "navbar-icon-links"],
    "show_prev_next": False,
    "show_version_warning_banner": True,
    # Footer
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version", "theme-version"],
    # Other
    "pygments_light_style": "default",
    "pygments_dark_style": "github-dark",
}
)

current_version = os.getenv("tag", "latest")

html_theme_options["switcher"] = dict(
    json_url="https://raw.githubusercontent.com/MartinPdeS/PyMieSim/documentation_page/version_switcher.json",
    version_match=current_version,
)

htmlhelp_basename = 'PyMieSimdoc'

latex_elements = {}


latex_documents = [
    (master_doc, 'PyMieSim.tex', 'PyMieSim Documentation',
     'Martin Poinsinet de Sivry-Houle', 'manual'),
]

man_pages = [
    (master_doc, 'pymiesim', 'PyMieSim Documentation',
     [author], 1)
]

texinfo_documents = [
    (master_doc, 'PyMieSim', 'PyMieSim Documentation',
     author, 'PyMieSim', 'One line description of project.',
     'Miscellaneous'),
]

epub_title = project

html_static_path = ['_static']
templates_path = ['_templates']
html_css_files = ['default.css']
epub_exclude_files = ['search.html']


# -- MyST --------------------------------------------------------------------
myst_enable_extensions = [
    # Enable fieldlist to allow for Field Lists like in rST (e.g., :orphan:)
    "fieldlist",
]
