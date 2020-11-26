
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Fields import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)


Detector = Photodiode(NA                = 0.2,
                      Source            = LightSource,
                      Npts              = 201,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)

Detector.PhiOffset = 45

Detector.Fourier.Plot('Polar')

Scat = Scatterer(Diameter      = 500e-9,
                 Source        = LightSource,
                 Index         = 1.4,
                 Meshes        = Detector.Meshes)

Scat.S1S2.Plot()

Scat.Stokes.Plot()

print(Detector.Coupling(Scatterer = Scat, Polarization='NoFiltered'))   # can be  all  -  Parallel  -  Perpendicular  -  Filtered  -  NoFiltered






# -
