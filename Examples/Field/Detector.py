
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Detector import Photodiode

LightSource = Source(Wavelength = 400e-9, Polarization = 0)

Detector = Photodiode(NA                = 0.8,
                      Source            = Source,
                      Sampling          = 1001,
                      GammaOffset       = 0,
                      PhiOffset         = 0)

Detector.Plot(num=100)





# -
