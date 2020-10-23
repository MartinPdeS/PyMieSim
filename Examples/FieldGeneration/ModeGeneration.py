
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""


from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode


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


LP01.PlotFields()

LP11.PlotFields()




# -
