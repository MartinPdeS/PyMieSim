
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.functions.converts import Angle2Direct, Direct2Angle
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Detector import Photodiode

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)

npts=101


LP11 = LPmode(Name        = 'LP11',
              Mode        = (1, 1),
              Source      = LightSource,
              Npts        = npts,
              NA          = 1)



Scat = Scatterer(Diameter    = 50e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Meshes      = LP11.FarField.Meshes)


print(LP11.Coupling(Scat, Mode='Centered'))

LP11.FarField.Plot()

Scat.FarField.Plot()

plt.show()
# -
