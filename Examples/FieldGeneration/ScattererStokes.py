
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter    = 200e-9,
                 Source      = LightSource,
                 Index       = 1.5,
                 Npts        = 201,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [-180, 180],
                 CacheTrunk  = None,
                 GPU         = False)


Scat.S1S2.Plot()


Scat.SPF.Plot()

Scat.Stokes.Plot()












# -
