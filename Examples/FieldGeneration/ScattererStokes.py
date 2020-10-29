
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

from PyMieCoupling.classes.Scattering import Scatterer


Scat = Scatterer(Diameter    = 200e-9,
                 Wavelength  = 400e-9,
                 Index       = 1.5,
                 Npts        = 201,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [-180, 180],
                 CacheTrunk  = None)


Scat.S1S2.Plot()

Scat.SPF.Plot()

Scat.Stokes.Plot()












# -
