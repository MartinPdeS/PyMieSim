# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
autodoc_mock_imports = ['numpy',
                        'pandas',
                        'typing',
                        'matplotlib',
                        'mayavi',
                        'scipy',
                        'tvtk',
                        'numpydoc',
                        'PyMieSim.LMT',
                        'PyMieSim.GLMT',
                        'PyMieSim._Coupling',
                        'PyMieSim._utils',
                        'beartype',
                        'PyMieSim.Fibonacci',
                        'sphinxcontrib.bibtex',
                        'sphinx.ext.autodoc']


PATH = os.path.join(os.path.dirname(__file__), '../..' )

sys.path.insert(0, PATH)


project   = 'PyMieSim'
copyright = '2021, Martin Poinsinet de Sivry-Houle'
author    = 'Martin Poinsinet de Sivry-Houle'


version    = ''

release    = '0.1.11'

extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    #'sphinxcontrib.bibtex',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.coverage',
    'sphinx.ext.autosectionlabel',
    'numpydoc'
]

numpydoc_show_class_members = False
# Add any paths that contain templates here, relative to this directory.

#bibtex_bibfiles = ['refs.bib']
#bibtex_encoding = 'latin'
#bibtex_default_style = 'unsrt'

templates_path = ['_templates']


# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

import os

#html_extra_path = ['../../tests/Examples/']
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
#pygments_style = None

pygments_style = 'sphinx'
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#




html_theme = 'sphinxdoc'








#html_theme_path = [ovs_sphinx_theme.get_theme_dir()]


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_static_path = ['_static']

#html_css_files = ['default.css']


html_logo = 'Logo.png'


# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'PyMieSimdoc'



html_css_files = [ 'default.css' ]

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'PyMieSim.tex', 'PyMieSim Documentation',
     'Martin Poinsinet de Sivry-Houle', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'pymiesim', 'PyMieSim Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'PyMieSim', 'PyMieSim Documentation',
     author, 'PyMieSim', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']

def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))

#list_files('../..')
# -- Extension configuration -------------------------------------------------
