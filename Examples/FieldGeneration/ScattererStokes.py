
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter    = 100e-9,
                 Source      = LightSource,
                 Index       = 1.5,
                 Npts        = 101,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [-180, 180])


Scat.S1S2.Plot()

Scat.SPF.Plot()

Scat.Stokes.Plot()








# -
