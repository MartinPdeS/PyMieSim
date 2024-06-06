#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import PyMieSim

__all__ = [
    'root_path',
    'project_path',
    'test_path',
    'static_doc_path',
    'lp_mode_path',
    'examples_path',
    'version_path',
    'validation_data_path',
    'doc_path',
    'logo_path',
    'doc_css_path'
]

root_path = Path(PyMieSim.__path__[0])

project_path = root_path.parents[0]

test_path = project_path.joinpath('tests')

static_doc_path = root_path.parents[0].joinpath('docs/images')

lp_mode_path = root_path.joinpath('lp_modes')

examples_path = root_path.joinpath('examples')

version_path = root_path.joinpath('VERSION')

validation_data_path = root_path.joinpath('validation_data')

doc_path = root_path.parents[0].joinpath('docs')

logo_path = doc_path.joinpath('images/logo.png')

doc_css_path = doc_path.joinpath('source/_static/default.css')

rtd_example = 'https://pymiesim.readthedocs.io/en/latest/gallery/index.html'

if __name__ == '__main__':
    for path_name in __all__:
        path = locals()[path_name]
        print(path)
        assert path.exists(), f"Path {path_name} do not exists"

# -
