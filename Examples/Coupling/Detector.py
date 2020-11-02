
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

npts = 201

GPU = True

Detector = Photodiode(NA                = 0.2,
                      Source            = LightSource,
                      Npts              = npts,
                      ThetaOffset       = 0,
                      PhiOffset         = 0,
                      GPU               = GPU)

Detector.PhiOffset = 45

Detector.Fourier.Plot('Polar')

Scat = Scatterer(Diameter      = 500e-9,
                 Source        = LightSource,
                 Index         = 1.4,
                 Meshes        = Detector.Meshes,
                 GPU           = GPU)

Scat.S1S2.Plot()

Scat.Stokes.Plot()








# -
