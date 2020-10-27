
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Photodiode

npts = 401

Detector = Photodiode(NumericalAperture = 0.3,
                      Wavelength        = 400e-9,
                      Npts              = npts,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)


Detector.PlotPolar()







# -
