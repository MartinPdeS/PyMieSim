#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from sphinx_gallery.sorting import FileNameSortKey

from PyMieSim.tools.directories import (
    logo_path,
    project_path,
    doc_css_path,
    version_path,
    examples_path
)

sys.path.insert(0, project_path)
sys.path.insert(0, project_path.joinpath('PyMieSim'))


def setup(app):
    app.add_css_file(str(doc_css_path))


autodoc_mock_imports = [
    'numpy',
    'matplotlib',
    'DataVisual',
    'MPSPlots',
    'scipy',
    'PyMieSim.Tools.measure',
    'PyMieSim.bin'
]

project = 'PyMieSim'
copyright = '2021, Martin Poinsinet de Sivry-Houle'
author = 'Martin Poinsinet de Sivry-Houle'
today_fmt = '%B %d, %Y'

with open(version_path, "r+") as f:
    version = release = f.read()


extensions = [
    'sphinx.ext.mathjax',
    'numpydoc',
    'pyvista.ext.plot_directive',
    'sphinx_gallery.gen_gallery',
]


try:
    import pyvista
    if sys.platform in ["linux", "linux2"]:
        pyvista.start_xvfb()  # Works only on linux system!
except ImportError:
    print('Could not load pyvista library for 3D renderin')


sphinx_gallery_conf = {
    'examples_dirs': [
        examples_path.joinpath('detector'),
        examples_path.joinpath('scatterer'),
        examples_path.joinpath('experiment/sphere'),
        examples_path.joinpath('experiment/cylinder'),
        examples_path.joinpath('experiment/coreshell'),
        examples_path.joinpath('validation')
    ],
    'gallery_dirs': [
        "Gallery/detector",
        "Gallery/scatterer",
        "Gallery/experiment/sphere",
        "Gallery/experiment/cylinder",
        "Gallery/experiment/coreshell",
        "Gallery/validation"
    ],
    'image_scrapers': ('matplotlib', 'pyvista'),
    'ignore_pattern': '/__',
    'plot_gallery': True,
    'thumbnail_size': [600, 600],
    'download_all_examples': False,
    'line_numbers': True,
    'remove_config_comments': True,
    'default_thumb_file': logo_path,
    'notebook_images': logo_path,
    'within_subsection_order': FileNameSortKey,
    'capture_repr': ('_repr_html_', '__repr__'),
    'nested_sections': True,
}

autodoc_default_options = {
    'members': False,
    'members-order': 'bysource',
    'undoc-members': False,
    'show-inheritance': True,
}

numpydoc_show_class_members = False

source_suffix = '.rst'

master_doc = 'index'

language = 'en'

exclude_patterns = []

pygments_style = 'monokai'

highlight_language = 'python3'

html_theme = 'sphinxdoc'

html_theme_options = {"sidebarwidth": 300}

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
