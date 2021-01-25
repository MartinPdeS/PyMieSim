
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""

from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Scattering import Scatterer, WMSample
from PyMieCoupling.classes.Detector import Photodiode, LPmode

LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)


LP11 = LPmode(Mode        = (0, 1),
              Sampling    = 401,
              NA          = 0.01,
              GammaOffset = 0,
              PhiOffset   = 0
              )



Scat = Scatterer(Diameter    = 1000e-9,
                 Source      = LightSource,
                 Index       = 1.5)

Sam = WMSample(g           = 0.8,
               lc          = 10*1e-7,
               D           = 2,
               Nc          = 10,
               Source      = LightSource)


Sam.Plot()
#LP11.Plot()

#Scat.Plot()




# -
