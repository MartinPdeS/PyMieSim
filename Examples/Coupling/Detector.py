
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode, LPmode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.utils import Source
import matplotlib.pyplot as plt

LightSource = Source(Wavelength   = 500e-9,
                     Polarization = 0)


Detector = Photodiode(Name              = 'Detector',
                      NA                = 0.5,
                      Source            = LightSource,
                      Samples           = 501,
                      GammaOffset       = 0,
                      PhiOffset         = 0)


LP11 = LPmode(Name        = 'LP11',
              Mode        = (1, 1),
              Source      = LightSource,
              Samples     = 501,
              NA          = 0.5,
              GammaOffset = 0,
              PhiOffset   = 0
              )


LP01 = LPmode(Name        = 'LP01',
              Mode        = (0, 1),
              Source      = LightSource,
              Samples     = 501,
              NA          = 0.5,
              GammaOffset = 0,
              PhiOffset   = 0
              )


Scat = Scatterer(Diameter      = 500e-9,
                 Source        = LightSource,
                 Index         = 1.4,
                 Meshes        = Detector.Meshes)

#Scat.S1S2.Plot()

#Scat.Parallel.Plot()

#Detector.Plot()

Coupling0 = Detector.Coupling(Scatterer = Scat, Mode='Centered')

Coupling1 = LP11.Coupling(Scatterer = Scat, Mode='Centered')

Coupling2 = LP01.Coupling(Scatterer = Scat, Mode='Centered')

print(Coupling0, Coupling1, Coupling2)

plt.show()

# -
