"""Shared dashboard defaults loaded from the bundled asset configuration."""

from __future__ import annotations

import json
from pathlib import Path


_DEFAULTS_PATH = Path(__file__).with_name("assets") / "default-settings.json"

with _DEFAULTS_PATH.open(encoding="utf-8") as _defaults_file:
    _defaults = json.load(_defaults_file)

DEFAULT_APPLICATION_SETTINGS = _defaults["application"]
DEFAULT_WORKSPACE_SETTINGS = _defaults["workspaces"]
DEFAULT_PLOT_SETTINGS = _defaults["plot"]
