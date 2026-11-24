"""Tests for the native extension installation preflight."""

from PyMieSim._native import REQUIRED_NATIVE_MODULES, missing_native_extensions


def test_all_required_native_extensions_are_discoverable():
    assert missing_native_extensions() == (), REQUIRED_NATIVE_MODULES
