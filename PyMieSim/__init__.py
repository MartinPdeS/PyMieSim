"""
PyMieSim package initialization.
"""

try:
    from ._version import version as __version__

except ImportError:
    __version__ = "0.0.0"

debug_mode = False

from PyMieSim.binary import interface_pint
from TypedUnit import ureg

interface_pint.set_ureg(ureg)


import PyMieSim.single.detector as _
import PyMieSim.single.scatterer as _
import PyMieSim.single.source as _
import PyMieSim.single.setup as _
import PyMieSim.mesh as _
import PyMieSim.coordinates as _
import PyMieSim.material as _
