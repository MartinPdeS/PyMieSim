
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 mode
_________________________________________________________
"""

from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling, MeanFieldCoupling


npts=200

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP01 = mode(fiber      = Fiber,
            LPmode     = (0, 1),
            wavelength = 400e-9,
            npts       = npts,
            PhiOffset = 10,
            ThetaOffset = 0,
            )

LP01.PlotFields()

Scat = Scatterer(diameter    = 200e-9,
                 wavelength  = 400e-9,
                 index       = 1.5,
                 npts        = 200,
                 Meshes = LP01.Meshes,
                 CacheTrunk = None)

Scat.PlotS1S2()

Scat.Field.PlotStokes(RectangleTheta=[-5,5],
                      RectanglePhi=[-5,5])














# -
