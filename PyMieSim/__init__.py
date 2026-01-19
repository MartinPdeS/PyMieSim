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