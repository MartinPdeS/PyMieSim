
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.classes.Scattering import FullScatterer, Scatterer
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.utils import Source

LightSource = Source(Wavelength   = 100e-9,
                     Polarization = 90)

Scat = Scatterer(Diameter    = 200e-9,
                 Source      = LightSource,
                 Acceptance  = 10,
                 Index       = 1.4,
                 Samples     = 2000,
                 PhiOffset   = 0,
                 GammaOffset = 30)

Scat.Plot()

Scat.S1S2.Plot()

Scat.SPF.Plot()

Scat.Meshes.Plot()

plt.show()





# -
