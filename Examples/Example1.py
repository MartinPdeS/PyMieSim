
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

Scat = Scatterer(diameter      = 500e-9,
                 wavelength    = 400e-9,
                 index         = 1.4,
                 npts          = npts,
                 Meshes = Detector.Meshes)

Scat.PlotS1S2()

Scat.Field.PlotStokes(RectangleTheta = [-5,5],
                      RectanglePhi   = [-5,5])

PointFieldCoupling(Detector = Detector,
                   Source   = Scat,
                   Field    = 'Parallel')










# -
