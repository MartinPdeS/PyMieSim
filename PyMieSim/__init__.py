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


from PyMieSim.single import detector as _
from PyMieSim.single import scatterer as _
from PyMieSim.single import source as _
from PyMieSim.single import setup as _
from PyMieSim.mesh import *
from PyMieSim.material import *