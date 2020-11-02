
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.Detector import Photodiode

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

npts = 401

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Detector = Photodiode(NumericalAperture = 0.3,
                      Source            = Source,
                      Npts              = npts,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)


Detector.PhiOffset = 20

Detector.Fourier.Plot('Polar')







# -
