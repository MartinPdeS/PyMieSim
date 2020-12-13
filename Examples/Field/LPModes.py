
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import fiber, LPmode
from PyMieCoupling.classes.Fields import Source
import matplotlib.pyplot as plt

npts=101

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 90)


LP11 = LPmode(Name        = 'LP11',
              Mode        = (1, 1),
              Source      = LightSource,
              Npts        = npts,
              ThetaOffset = 0,
              PhiOffset   = 0,
              NA          = 0.4)

LP01 = LPmode(Name          = 'LP01',
              Mode          = (0, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              NA            = 0.4)


print(LP01.FarField.Spherical.shape)

LP01.PhiOffset = 10

LP01.ThetaOffset = 60

LP01.NearField.Plot()

LP01.FarField.Plot()

plt.show()
# -
