
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""


from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.Scattering import Scatterer


LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

npts=21

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = LPmode(Fiber         = Fiber,
              Mode          = (1, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Magnification = 1)

LP01 = LPmode(Fiber         = Fiber,
              Mode          = (0, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Magnification = 1)

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter    = 500e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Meshes      = LP01.Meshes)


print( LP01.Coupling(Scat) )

print( LP11.Coupling(Scat) )



# -
