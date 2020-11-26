
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import fiber, LPmode
from PyMieCoupling.classes.Fields import Source

npts=101

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)


Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = LPmode(Fiber       = Fiber,
              Name        = 'LP11',
              Mode        = (1, 1),
              Source      = LightSource,
              Npts        = npts,
              ThetaOffset = 0,
              PhiOffset   = 0,
              NA          = 0.4)

LP01 = LPmode(Fiber         = Fiber,
              Name          = 'LP01',
              Mode          = (0, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              NA            = 0.4)


LP01.Field.Plot()

LP11.Field.Plot()

LP01.Fourier.Plot()

LP11.Fourier.Plot()

LP01.Fourier.Plot('Polar')

LP11.Fourier.Plot('Polar')


# -
