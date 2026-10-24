#!/usr/bin/env python3
"""Delegate PyMieSim release operations to martinpdes-maintenance."""

import sys
from pathlib import Path

from martinpdes_maintenance.cli import main

if __name__ == "__main__":
    raise SystemExit(main(["--root", str(Path(__file__).resolve().parents[1]), "check", *sys.argv[1:]]))
