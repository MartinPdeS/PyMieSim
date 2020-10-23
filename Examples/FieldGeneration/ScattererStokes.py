
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling, MeanFieldCoupling


npts=200

Scat = Scatterer(diameter    = 200e-9,
                 wavelength  = 400e-9,
                 index       = 1.5,
                 npts        = 201,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [-180, 180],
                 CacheTrunk  = None)

Scat.PlotS1S2()

Scat.Field.PlotStokes(RectangleTheta = [-5,5],
                      RectanglePhi   = [-5,5])














# -
