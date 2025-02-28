from PyMieSim.units import Quantity  # God knows why removing this line make the test fails on mac x86_64! DO NOT REMOVE!

try:
    from ._version import version as __version__

except ImportError:
    __version__ = "0.0.0"

debug_mode = False

# -
