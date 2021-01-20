
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.functions.converts import Angle2Direct, Direct2Angle
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Scattering import Scatterer, Sample
from PyMieCoupling.classes.Detector import Photodiode

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 90)


LP11 = LPmode(Name        = 'LP11',
              Mode        = (0, 1),
              Source      = LightSource,
              Sampling    = 401,
              NA          = 0.01,
              GammaOffset = 0,
              PhiOffset   = 0
              )



Scat = Scatterer(Diameter    = 1000e-9,
                 Source      = LightSource,
                 Index       = 1.5,
                 Meshes      = LP11.Meshes)

Sam = Sample(g           = 0.8,
             lc          = 10*1e-7,
             D           = 2,
             Nc          = 10,
             Source      = LightSource,
             Meshes      = LP11.Meshes)


Sam.Plot()
#LP11.Plot()

#Scat.Plot()

plt.show()
# -
