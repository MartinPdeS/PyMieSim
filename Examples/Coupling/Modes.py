
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""



import numpy as np
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling


npts=201

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = mode(fiber       = Fiber,
            LPmode      = (1, 1),
            wavelength  = 400e-9,
            npts        = npts,
            ThetaOffset = 0,
            PhiOffset   = 0
            )

LP01 = mode(fiber       = Fiber,
            LPmode      = (0, 1),
            wavelength  = 400e-9,
            npts        = npts,
            ThetaOffset = 0,
            PhiOffset   = 0
            )



LP01.magnificate(magnification=2.)

LP11.magnificate(magnification=2.)

Scat = Scatterer(diameter    = 500e-9,
                 wavelength  = 400e-9,
                 index       = 1.4,
                 Meshes      = LP01.Meshes
                 )

LP01Perp, LP01Para = PointFieldCoupling(Detector = LP01,
                                        Source   = Scat,)

LP11Perp, LP11Para = PointFieldCoupling(Detector = LP11,
                                        Source   = Scat,)








# -
