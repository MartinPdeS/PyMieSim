
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Detector import Photodiode

LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

Detector = Photodiode(NA                = 0.8,
                      Sampling          = 1001,
                      GammaOffset       = 0,
                      PhiOffset         = 0)

Detector.Plot()





# -
