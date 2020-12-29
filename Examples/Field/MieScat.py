
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.classes.Scattering import FullScatterer, Scatterer
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.utils import Source

LightSource = Source(Wavelength   = 950e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter    = 5000e-9,
                 Source      = LightSource,
                 Acceptance  = 20,
                 Index       = 1.4,
                 Samples     = 500,
                 PhiOffset   = 20,
                 GammaOffset = 60)

fig0, fig1 = Scat.Plot(scatter=False)

Scat.S1S2.Plot()

Scat.SPF.Plot()

Scat.Meshes.Plot()
#fig0.savefig("lambert_5000.png", dpi=300)

plt.show()





# -
