
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling


npts = 201

Detector = Photodiode(NumericalAperture = 0.2,
                      Wavelength        = 400e-9,
                      Npts              = npts,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)

Detector.PlotPolar()

Scat = Scatterer(Diameter      = 500e-9,
                 Wavelength    = 400e-9,
                 Index         = 1.4,
                 Meshes        = Detector.Meshes)

Scat.PlotS1S2()

Scat.Field.PlotStokes(RectangleTheta = [-5,5],
                      RectanglePhi   = [-5,5])

Perp, Para = PointFieldCoupling(Detector = Detector,
                                Source   = Scat)







# -
