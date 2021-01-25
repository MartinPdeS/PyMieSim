
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.utils import Source


LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)


LP11 = LPmode(Mode         = (1, 1,'h'),
              Sampling     = 201,
              NA           = 0.8,
              GammaOffset  = 0,
              PhiOffset    = 0,
              CouplingMode = 'Centered')


#LP11.NA = 1

#LP11.Meshes.Plot()
LP11.Plot()

# --
