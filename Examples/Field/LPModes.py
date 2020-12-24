
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.utils import Source
import matplotlib.pyplot as plt

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 90)


LP11 = LPmode(Name        = 'LP11',
              Mode        = (1, 1),
              Source      = LightSource,
              Samples     = 201,
              NA          = 0.1,
              GammaOffset = 0,
              PhiOffset   = 0
              )

LP11.NA = 0.8

LP11.Meshes.Plot()

plt.show()
# --
