
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""

from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Misc import Source

npts=201

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = LPmode(Fiber       = Fiber,
              Name        = 'LP11',
              Mode        = (1, 1),
              Wavelength  = 400e-9,
              Npts        = npts,
              ThetaOffset = 0,
              PhiOffset   = 0
            )

LP01 = LPmode(Fiber       = Fiber,
              Name        = 'LP01',
              Mode      = (0, 1),
              Wavelength  = 400e-9,
              Npts        = npts,
              ThetaOffset = 0,
              PhiOffset   = 0
            )


LP01.magnificate(Magnification=2.)

LP11.magnificate(Magnification=2.)

LP01.PlotPolar()

LP11.PlotPolar()

LP01.PlotDirectSpace()




# -
