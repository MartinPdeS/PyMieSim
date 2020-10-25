
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import Detector
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling


npts = 201

Detector = Detector(size       = 50e-6,
                    wavelength = 400e-9,
                    shape      = 'circle',
                    npts       = npts,
                    ThetaOffset = 0,
                    PhiOffset = 0)

Detector.magnificate(magnification=1.5)

Detector.PlotFields()











# -
