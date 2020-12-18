
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.classes.Scattering import FullScatterer
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.utils import Source

LightSource = Source(Wavelength   = 100e-9,
                     Polarization = 0)

Scat = FullScatterer(Diameter    = 10e-9,
                     Source      = LightSource,
                     Index       = 1.4,)




Scat.FarField.Plot()

Scat.S1S2.Plot()

Scat.SPF.Plot()

Scat.Meshes.Plot()

plt.show()





# -
