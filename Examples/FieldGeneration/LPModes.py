
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""

from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Misc import Source

npts=101

cuda = True

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
              cuda        = cuda)

LP01 = LPmode(Fiber       = Fiber,
              Name        = 'LP01',
              Mode        = (0, 1),
              Source      = LightSource,
              Npts        = npts,
              ThetaOffset = 0,
              PhiOffset   = 0,
              cuda        = cuda)



LP01.Magnificate(Magnification=2.0)

#LP11.Magnificate(Magnification=2.0)



LP01.Field.Plot('Real')

LP11.Field.Plot('Imag')

LP01.Fourier.Plot('Polar')

LP11.Fourier.Plot('Polar')




# -
