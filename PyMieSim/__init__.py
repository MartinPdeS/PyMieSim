try:
    from ._version import version as __version__

except ImportError:
    __version__ = "0.0.0"

debug_mode = False

import PyMieSim.units as _
import PyMieSim.material as _

