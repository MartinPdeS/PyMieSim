
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.utils import Source

LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

Scat = Scatterer(Diameter    = 500e-9,
                 Source      = LightSource,
                 Index       = 1.4)

#Field = Scat.Field(Num=150)
SPF = Scat.SPF(Num=80)
#S1S2 = Scat.S1S2(Num=80)

SPF.Plot()




# -
