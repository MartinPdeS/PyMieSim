
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Fields import Source
import matplotlib.pyplot as plt

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 90)


Detector = Photodiode(NA                = 0.2,
                      Source            = LightSource,
                      Npts              = 201,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)




Scat = Scatterer(Diameter      = 500e-9,
                 Source        = LightSource,
                 Index         = 1.4,
                 Meshes        = Detector.FarField.Meshes)

#Scat.S1S2.Plot()

#Scat.Stokes.Plot()

Scat.Field.Plot()

#Detector.FarField.Plot()

print(Detector.Coupling(Scatterer = Scat, Mode='Centered'))

plt.show()

# -
