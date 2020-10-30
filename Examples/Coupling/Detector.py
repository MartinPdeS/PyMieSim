
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.functions.couplings import PointFieldCoupling

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

npts = 201

Detector = Photodiode(NumericalAperture = 0.2,
                      Source            = LightSource,
                      Npts              = npts,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)

Detector.PlotPolar()

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter      = 500e-9,
                 Source        = LightSource,
                 Index         = 1.4,
                 Meshes        = Detector.Meshes)

Scat.S1S2.Plot()

Scat.Stokes.Plot()

Perp, Para = PointFieldCoupling(Detector = Detector,
                                Source   = Scat)







# -
