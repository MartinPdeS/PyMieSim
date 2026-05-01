"""Top-level package for PyMieSim."""

try:
    from ._version import version as __version__

except ImportError:
    __version__ = "0.0.0"

import PyMieSim.units as _
import PyMieSim.material as _
import PyMieSim.polarization as _
