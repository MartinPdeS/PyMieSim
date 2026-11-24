"""Tests for the native extension installation preflight."""

import runpy
import PyMieSim
from PyMieSim import directories
from PyMieSim._native import REQUIRED_NATIVE_MODULES, missing_native_extensions


def test_all_required_native_extensions_are_discoverable():
    assert missing_native_extensions() == (), REQUIRED_NATIVE_MODULES


def test_validation_data_with_split_editable_package(monkeypatch, tmp_path):
    """Reference data lives beside the Python sources in a split editable install."""
    monkeypatch.setattr(PyMieSim, "__path__", [str(tmp_path), *PyMieSim.__path__])
    paths = runpy.run_path(directories.__file__)
    assert (paths["validation_data_path"] / "pymiescatt" / "validation_sphere.csv").is_file()
