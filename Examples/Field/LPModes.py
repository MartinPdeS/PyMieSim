
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.utils import Source
import matplotlib.pyplot as plt

LightSource = Source(Wavelength = 400e-9, Polarization = 90)


LP11 = LPmode(Name        = 'LP11',
              Mode        = (1, 1),
              Source      = LightSource,
              Samples     = 501,
              NA          = 0.8,
              GammaOffset = 250,
              PhiOffset   = 90
              )

#LP11.NA = 1

#LP11.Meshes.Plot()
LP11.Plot()
plt.show()
# --
