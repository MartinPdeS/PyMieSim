
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.functions.converts import Angle2Direct, Direct2Angle
from PyMieCoupling.classes.Detector import fiber, LPmode
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Scattering import Scatterer

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 30)

npts=201


LP11 = LPmode(Mode          = (1, 1),
              Orientation   = 'v',
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 10,
              Filter        = 0,
              NA            = 0.8)


LP01 = LPmode(Mode          = (0, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Filter        = 0,
              NA            = 0.8,
              Orientation   = 'v')


Scat = Scatterer(Diameter    = 15000e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Meshes      = LP11.FarField.Meshes)


print(LP11.Coupling(Scat, Mode='Centered'))

LP11.FarField.Plot()

Scat.Field.Plot()

plt.show()
# -
