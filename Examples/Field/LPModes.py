
"""
_________________________________________________________
Plottings of the far-field for LP modes.
_________________________________________________________
"""

from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.utils import Source
import matplotlib.pyplot as plt

npts=101

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 90)


LP11 = LPmode(Name        = 'LP11',
              Mode        = (1, 1),
              Source      = LightSource,
              Npts        = npts,
              NA          = 1)

LP01 = LPmode(Name          = 'LP01',
              Mode          = (0, 1),
              Source        = LightSource,
              Npts          = npts,
              NA            = 1)


#LP01.PhiOffset = 10

#LP01.ThetaOffset = 60

#LP01.NearField.Plot()

LP01.FarField.Plot()

plt.show()
# -
