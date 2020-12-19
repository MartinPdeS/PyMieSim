
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.utils import Source
import matplotlib.pyplot as plt

LightSource = Source(Wavelength   = 1000e-9,
                     Polarization = 90)


Detector = Photodiode(NA                = 0.8,
                      Source            = LightSource,
                      Npts              = 201,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)




Scat = Scatterer(Diameter      = 50e-9,
                 Source        = LightSource,
                 Index         = 1.4,
                 Meshes        = Detector.FarField.Meshes)

Scat.S1S2.Plot()

Scat.FarField.Plot()

Detector.FarField.Plot()

print(Detector.Coupling(Scatterer = Scat, Mode='Centered'))

plt.show()

# -
