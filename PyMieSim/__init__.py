"""PyMieSim package initialization.

This module exposes the package version and sets up unit handling.  Importing
``Quantity`` here ensures proper initialization on all platforms.  Removing this
import causes unit tests to fail on macOS systems.
"""

from PyMieSim.units import Quantity  # required for unit registration on macOS

try:
    from ._version import version as __version__

except ImportError:
    __version__ = "0.0.0"

debug_mode = False
