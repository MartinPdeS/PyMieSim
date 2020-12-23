
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.utils import Source
import matplotlib.pyplot as plt

LightSource = Source(Wavelength   = 100e-9,
                     Polarization = 0)


Detector = Photodiode(NA                = 0.9,
                      Source            = Source,
                      Samples           = 1001,
                      GammaOffset       = 0,
                      PhiOffset         = 0)



Scat = Scatterer(Diameter      = 500e-9,
                 Source        = LightSource,
                 Index         = 1.4,
                 Meshes        = Detector.Meshes)

#Scat.S1S2.Plot()

#Scat.Parallel.Plot()

#Detector.Plot()

print(Detector.Coupling(Scatterer = Scat, Mode='Centered'))

plt.show()

# -
