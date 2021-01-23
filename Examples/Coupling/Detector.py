
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode, LPmode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.utils import Source


LightSource = Source(Wavelength = 500e-9, Polarization = 0)


Detector = Photodiode(NA                = 0.5,
                      Source            = LightSource,
                      Sampling          = 501,
                      GammaOffset       = 0,
                      PhiOffset         = 0,
                      CouplingMode      = 'Centered')

Scat = Scatterer(Diameter      = 1000e-9,
                 Source        = LightSource,
                 Index         = 1.4,
                 Meshes        = Detector.Meshes)


Detector.Plot()

Coupling0 = Detector.Coupling(Scatterer = Scat)


# -
