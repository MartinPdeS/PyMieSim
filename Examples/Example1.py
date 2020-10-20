
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from src.classes.Detector import Detector
from src.classes.Scattering import Scatterer
from src.functions.couplings import PointFieldCoupling



npts = 201

Detector = Detector(size       = 50e-6,
                    wavelength = 400e-9,
                    shape      = 'circle',
                    npts       = npts)

#Detector.magnificate(magnification=1.5)

Detector.PlotFields()

Scat = Scatterer(diameter      = 500e-9,
                 wavelength    = 400e-9,
                 index         = 1.4,
                 npts          = npts,
                 ThetaBound    = [-180,180],
                 ThetaOffset   = 0,
                 PhiBound      = [-180,180],
                 PhiOffset     = 0)

Scat.PlotS1S2()

Scat.Field.PlotStokes(RectangleTheta = [-5,5],
                      RectanglePhi   = [-5,5])

PointFieldCoupling(Detector = Detector,
                   Source   = Scat.Field.Parallel,
                   Mesh     = Scat.Meshes)










# -
