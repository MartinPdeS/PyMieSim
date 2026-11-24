#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path

__all__ = ["root_path", "validation_data_path", "doc_path", "logo_path", "doc_css_path"]

# Editable builds can put native extensions first in the package search path.
root_path = Path(__file__).resolve().parent

validation_data_path = root_path.joinpath("validation_data")

doc_path = root_path.parents[0].joinpath("docs")

logo_path = doc_path.joinpath("images/logo.png")

doc_css_path = doc_path.joinpath("source/_static/default.css")

if __name__ == "__main__":
    for path_name in __all__:
        path = locals()[path_name]
        print(path)
        assert path.exists(), f"Path {path_name} do not exists"

# -
