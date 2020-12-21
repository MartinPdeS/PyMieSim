
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
                     Polarization = 0)
"""
Scat = FullScatterer(Diameter    = 100e-9,
                     Source      = LightSource,
                     Index       = 1.4)
"""
Scat = Scatterer(Diameter    = 100e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Samples     = 1000,
                 PhiOffset   = 40)

Scat.Parallel.Plot()

Scat.S1S2.Plot()

Scat.SPF.Plot()

Scat.Meshes.Plot()

plt.show()





# -
