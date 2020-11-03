
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.Detector import Photodiode

npts = 401

cuda = True

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Detector = Photodiode(NA                = 0.3,
                      Source            = Source,
                      Npts              = npts,
                      ThetaOffset       = 0,
                      PhiOffset         = 0,
                      cuda              = cuda)


Detector.PhiOffset = 20

Detector.Fourier.Plot('Real')

Detector.Fourier.Plot('Imag')

Detector.Fourier.Plot('Polar')







# -
